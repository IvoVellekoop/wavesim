classdef(Abstract) WaveSimBase < Simulation
    % Base class for the simulation of 2-D or 3-D wave equations using
    % the modified Born series approach
    % This class is overridden by WaveSim (for scalar simulations)
    % and WaveSimVector (for vector simulations).
    %
    % The actual simulation is performed by calling the functions
    % 'mix' and 'propagate'. Only the 'propagate' function is different
    % for vector or scalar waves.
    %
    % Ivo M. Vellekoop 2014-2018
    
    properties
        gamma;     % potential array used in simulation
        epsilon;   % convergence parameter
        k02e;      % precomputed constants (pre-divided by sqrt epsilon)
        mix;       % function handle to function performing the mixing step
        %% wiggle properties (might be moved to separate class?)
        wiggles;
        % cell array containing all the instructions for the
        % anti-wraparound and anti-aliasing algorithm (the wiggle algorithms)
        boundary_wiggle = true;         
        % 'true' indicates that the anti-wraparound algorithm is used on
        % all edges with non-zero boundary width. Note: zero-width
        % boundaries are treated as periodic boundaries, and the
        % anti-wraparound algorithm is disabled by default for these
        % boundaries. To override this default behavior, you can explicitly
        % pass a 3x1 logical vector (e.g. [true false false]) to indicate
        % which boundaries to 'wiggle'
        medium_wiggle = false; 
        % true indicates that the anti-aliasing algorithm is used in all 
        % dimensions. Algorithm can also exclusively enabled in a single or 
        % two dimensions by passing a 3x1 logical vector. 
        Nwiggles; 
        % number of wiggles performed (all combination of boundary_wiggle 
        % and medium_wiggle)
        
        %% diagnostics and feedback
        epsilonmin = 3; % minimum value to avoid divergence when simulating empty medium   
    end
    methods(Abstract)
        propagate(obj); % function performing the propagation step        
    end
    methods
        function obj = WaveSimBase(sample, options)
            %% Constructs a wave simulation object
            %	sample = SampleMedium object
            %   options.lambda = free space wavelength (same unit as pixel_size, e. g. um)
            %   options.epsilon = convergence parameter (leave empty unless forcing a specific value)            
            options.boundary_wiggle = options.wiggle; % temporarily fix for change in syntax
            obj@Simulation(sample, options);
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use
            
            %% determine wavenumbers
            % the optimal k_0 follows is given by n_center (see Medium)
            k0 = 2*pi/obj.lambda;                   % wavenumber in empty medium
            k0c = sqrt(sample.e_r_center) * k0; % wavenumber for center value refractive index
            
            %% Determine epsilon and gamma    
            % determine optimal epsilon (if value is not given in options) 
            if ~isfield(options, 'epsilon') % setting epsilon in options forces a specific value, may not converge
                obj.epsilon = obj.calculate_epsilon(sample.e_r,k0,k0c);
            end
            obj.iterations_per_cycle = obj.lambda /(2*k0c/obj.epsilon); %divide wavelength by pseudo-propagation length    
            
            % calculate gamma for all submedia
            obj.gamma = cell(sample.n_media,1);
            for i_medium = 1:sample.n_media                
                V = sample.e_r{i_medium}*k0^2-k0c^2 - 1.0i*obj.epsilon; % calculate potential map                   
                V = Medium.apply_filters(V,sample.filters); % apply sample filters to edge of potential map
                obj.gamma{i_medium} = obj.data_array(1.0i / obj.epsilon * V); % calculate gamma and convert to desired data type
            end
            obj.gamma = obj.gamma{1};           % medium_wiggle currently still disabled (still under development)
            
            %% convert to single or double precision, and put on gpu if
            % needed. 
            obj.k02e = obj.data_array(k0c^2/obj.epsilon + 1.0i);
            obj.epsilon = obj.data_array(obj.epsilon);      

            %% define function for the mixing step:
            % Ediff = (1-gamma) Ediff + gamma^2 [G Ediff]
            % (where [G Ediff] is the propagated field)
            %
            mix = @(Ediff, Eprop, gamma) (1-gamma) .* Ediff - 1.0i * gamma.^2 .* Eprop;
            
            if obj.gpu_enabled
                % For performance reasons, we use arrayfun when using the
                % GPU. This way we prevent unnessecary memory allocations
                % and memory access. Unfortunately, using arrayfun on the
                % CPU is extremely inefficient (tested in Matlab 2017a)
                %
                obj.mix = @(Eold, Enew, gamma) arrayfun(mix, Eold, Enew, gamma);
            else
                obj.mix = mix;
            end
            
            %% calculate wiggle descriptors
            [obj.wiggles, obj.boundary_wiggle, obj.medium_wiggle, obj.Nwiggles] = obj.compute_wiggles();    
        end
        
        function state = run_algorithm(obj, state)
            % paper: 
            %   E = [1+M+M^2+...] \gamma G S
            % iterate:
            %   E_{0}   = 0
            %   E_{k+1} = M E_k + \gamma G S
            %   M = \gamma G V - \gamma + 1
            %
            % Now instead accumulate all in one buffer. dE_{k} are the terms in the Born series:
            %   dE_{k} = M^k \gamma G S
            % 
            % or, recursively:
            %   dE_{1} = \gamma G S
            %   dE_{k} = M dE_{k-1}
            %
            % and we accumulate all in one buffer to find E
            %   E = dE_{1} + dE_{2} + ...
            %
            % we now substitute V = \epsilon/i \gamma and define
            % G' = \epsilon G
            % so that M = -i \gamma G' \gamma - \gamma + 1
            %
            %   dE_{1} = \gamma G' i/\epsilon S
            %   dE_{k} = (-i \gamma G' \gamma - \gamma + 1) dE_{k-1}
            %
            % A further optimization can be done by replacing dE = dE' / \gamma
            %
            %   dE'_{1} = \gamma^2 G' i/\epsilon S
            %   dE'_{k} = (-i \gamma^2 G' - \gamma + 1) dE'_{k-1}
            %           = -i \gamma^2 G' dE'_{k-1}  + (1-\gamma) dE'_{k-1}
            %
            % Which simplifies to
            %
            %   dE'_{0} = 0
            %   dE'_{1} = -i\gamma^2 G' [dE'_{k-1} + i/\epsilon S]  + (1-\gamma) dE'_{k-1}
            %   dE'_{k} = -i\gamma^2 G'  dE'_{k-1}                  + (1-\gamma) dE'_{k-1}
            
            %   These iterations are implemented as:
            %   1. a propagation step:   Eprop   =  G' [E_diff + i/\epsilon S] {only add S in 1st iteration}
            %   2. a mixing step:        E_diff  => (1-\gamma) E_diff - i \gamma^2 Eprop  
            %   3. an accumulation step: E       => E + E_diff
            %   
            %   And, after all iterations:
            %       E       => E / \gamma
            
            %% Allocate memory for calculations
            state.E = obj.data_array([], obj.N);
            state.Ediff = obj.data_array([], obj.N);
            
            %% simulation iterations
            while state.has_next
                wiggle = obj.wiggles{mod(state.it, obj.Nwiggles) + 1}; 
                if state.it <= obj.Nwiggles % During the first few iterations: add source term
                    Etmp = state.source.add_to(state.Ediff, 1.0i / obj.epsilon / obj.Nwiggles);
                else
                    Etmp = state.Ediff;
                end
                
                Etmp = obj.propagate(Etmp, wiggle);
                state.Ediff = obj.mix(state.Ediff, Etmp, obj.gamma);
                state.E = state.E + state.Ediff;
                
                if state.calculate_energy
                    state.last_step_energy = Simulation.energy( obj.crop_field(state.Ediff) );                  
                end
                
                can_terminate = mod(state.it, obj.Nwiggles) == 0; %only stop after multiple of Nwiggles iterations
                state = next(obj, state, can_terminate);                
            end
            
            % divide field by gamma to convert E' -> E and crop field to 
            % remove boundary layers 
            state.E = state.E ./ obj.gamma;
            state.E = obj.crop_field(state.E);
        end
    end
    methods(Access=private)
        function epsilon = calculate_epsilon(obj, e_r, k0, k0c)
            % function used to calculate the optimal convergence parameter 
            % epsilon used by wavesim based on maximum real value of 
            % potential map   
            
            % determine largest absolute potential fluctuation of all submedia 
            % to guarantee convergence
            V_tot = cell2mat(e_r)*k0^2-k0c^2;
            [Vabs_max, max_index] = max(abs(V_tot(:)));
            
            % if Vabs_max corresponds to a point with only absorption
            % and we use epsilon = Vabs_max, then at that point V will be
            % 0, and the method does not work. In this special case, add
            % a small offset to epsilon
            if real(V_tot(max_index)) < 0.05 * Vabs_max
                Vabs_max = Vabs_max * 1.05;
            end
            epsilon = max(Vabs_max, obj.epsilonmin); % minimum value epsilonmin to avoid divergence when simulating empty medium
        end
                
        function Ecrop = crop_field(obj,E)
            % Removes the boundary layer from the simulated field by
            % cropping field dataset
            Ecrop = E(obj.output_roi(1,1):obj.output_roi(2,1),...
                      obj.output_roi(1,2):obj.output_roi(2,2),...
                      obj.output_roi(1,3):obj.output_roi(2,3),...
                      obj.output_roi(1,4):obj.output_roi(2,4));
        end
        
        %%% Wiggle methods (move to separate class?)
        function [wiggle_descriptors,boundary_wiggle,medium_wiggle,Nwiggles] = compute_wiggles(obj) 
            % Decides which borders to wiggle and returns wiggled coordinates
            % for those borders

            % initialize boundary wiggle
            % true -> wiggle all non-periodic boundaries
            % [true, false, true] -> wiggle only 1st and 3rd dimension
            % false -> same as [false, false, false]
            boundary_wiggle = obj.boundary_wiggle(:);
            if numel(boundary_wiggle) == 1
                if boundary_wiggle == true
                    boundary_wiggle = ~obj.grid.periodic(:);
                else
                    boundary_wiggle = false(3,1);
                end
            end
            
            % initialize medium_wiggle vector
            medium_wiggle = obj.medium_wiggle(:);
            if numel(medium_wiggle) == 1 % either true or false in all dimensions
                medium_wiggle = logical(medium_wiggle*ones(3,1));
            end
                            
            % generate matrix with all posible wiggle combinations (combine 
            % boundary_wiggle and medium_wiggle)
            wiggles_3D = [-1,  1, -1, -1,  1,  1, -1, 1;...
                          -1, -1,  1, -1,  1, -1,  1, 1;...
                          -1, -1, -1,  1, -1,  1,  1, 1]; % all possible permutations in 3D
            all_wiggles = zeros(6,2^6); % [y;x;z;ky;kx;kz] in all 64 different permutations 
            for w = 1:8
                all_wiggles(:,(w-1)*8+(1:8)) = [wiggles_3D; repmat(wiggles_3D(:,w),1,8)];
            end
             
            % determine wiggle directions
            wiggles_combined = [boundary_wiggle;medium_wiggle]';
            [~,idx] = unique(all_wiggles' .* wiggles_combined, 'rows');
            boundary_wiggles = all_wiggles(1:3,sort(idx)).* boundary_wiggle;
            medium_wiggles = all_wiggles(4:6,sort(idx)).* medium_wiggle;
            
            % calculate phase ramps and coordinates for every
            % boundary_wiggle and medium_wiggle
            Nwiggles = numel(idx);
            wiggle_descriptors = cell(Nwiggles,1); % pre-allocate memory
            for w_i=1:Nwiggles 
                wiggle_descriptors{w_i} = obj.wiggle_descriptor(boundary_wiggles(:, w_i), medium_wiggles(:, w_i));
            end
        end
        
        function wd = wiggle_descriptor(obj, b_wigg, m_wigg)
            % Constructs shifted coordinates and phase ramps for a wiggle
            % steps. (b_wigg: boundary_wiggle, m_wigg: medium_wiggle)
            wd = struct;
            % construct coordinates, shift quarter of a pixel when wiggling
            % (required for calculating wiggled Green's function). Pre-scale 
            % Fourier coordinates to optimize the propagation functions a bit
            wd.pxe = obj.data_array((obj.grid.px_range - obj.grid.dpx * b_wigg(2)/4)/sqrt(obj.epsilon));
            wd.pye = obj.data_array((obj.grid.py_range - obj.grid.dpy * b_wigg(1)/4)/sqrt(obj.epsilon));
            wd.pze = obj.data_array((obj.grid.pz_range - obj.grid.dpz * b_wigg(3)/4)/sqrt(obj.epsilon));
            
            % construct real space phase gradients to compensate for the
            % pixel shift in k_space (required for boundary_wiggle)
            wd.gx = obj.data_array(exp(1.0i * pi/2 * b_wigg(2) * obj.grid.x_range / obj.grid.dx / length(obj.grid.x_range)));
            wd.gy = obj.data_array(exp(1.0i * pi/2 * b_wigg(1) * obj.grid.y_range / obj.grid.dx / length(obj.grid.y_range)));
            wd.gz = obj.data_array(exp(1.0i * pi/2 * b_wigg(3) * obj.grid.z_range / obj.grid.dx / length(obj.grid.z_range)));
            
            % construct k-space phase gradients to compensate for the
            % pixel shift in real space (required for medium_wiggle)
            wd.gpx = obj.data_array(exp(1.0i * pi/2 * m_wigg(2) * obj.grid.px_range / obj.grid.dpx / length(obj.grid.px_range)));
            wd.gpy = obj.data_array(exp(1.0i * pi/2 * m_wigg(1) * obj.grid.py_range / obj.grid.dpy / length(obj.grid.py_range)));
            wd.gpz = obj.data_array(exp(1.0i * pi/2 * m_wigg(3) * obj.grid.pz_range / obj.grid.dpz / length(obj.grid.py_range)));      
        end
    end
end


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
        no_wiggle; % holds coordinates for non-wiggled case
        wiggle = true;         
        % 'true' indicates that the anti-wraparound algorithm is used on
        % all edges with non-zero boundary width. Note: zero-width
        % boundaries are treated as periodic boundaries, and the
        % anti-wraparound algorithm is disabled by default for these
        % boundaries. To override this default behavior, you can explicitly
        % pass a 1x3 logical vector (e.g. [true false false]) to indicate
        % which boundaries to 'wiggle'
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
            obj.gamma = obj.gamma{1};
            
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
            
            %% calculate wiggle coefficients
            [obj.wiggle, obj.no_wiggle] = obj.compute_wiggles(obj.wiggle);    
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
            Nwiggle = size(obj.wiggle, 2);
            while state.has_next
                wigg = obj.wiggle(mod(state.it, Nwiggle) + 1); 
                if state.it <= Nwiggle % During the first few iterations: add source term
                    Etmp = state.source.add_to(state.Ediff, 1.0i / obj.epsilon / Nwiggle);
                else
                    Etmp = state.Ediff;
                end
                Etmp = obj.propagate(Etmp, wigg);
                state.Ediff = obj.mix(state.Ediff, Etmp, obj.gamma);
                state.E = state.E + state.Ediff;
                
                if state.calculate_energy
                    state.last_step_energy = Simulation.energy( obj.crop_field(state.Ediff) );                  
                end
                
                can_terminate = mod(state.it, Nwiggle) == 0; %only stop after multiple of Nwiggle iterations
                state = next(obj, state, can_terminate);                
            end
            
            state.E = state.E ./ obj.gamma;
            state.rel_error = obj.calculate_rel_error(state);
            
            % crop field to remove boundary layers
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
        
        function [wiggle_descriptors, no_wiggle] = compute_wiggles(obj, wiggle_option) 
            %% Decides which borders to wiggle and returns wiggled coordinates
            % for those borders

            % calculate descriptor without wiggling
            no_wiggle = obj.wiggle_descriptor([0;0;0]);

            % decide which borders to wiggle
            % true -> wiggle all non-periodic boundaries
            % [true, false, true] -> wiggle only 1st and 3rd dimension
            % false -> same as [false, false, false]
            %
            if numel(wiggle_option) == 1 && wiggle_option == true %'auto'
                wiggle_option = ~obj.grid.periodic; 
            end
            
            if ~any(wiggle_option)
                wiggle_descriptors = no_wiggle;
                return;
            end
            
            %determine wiggle directions
            wiggles = ...
                [0, 1, 0, 0, 1, 1, 0, 1;...
                 0, 0, 1, 0, 1, 0, 1, 1;...
                 0, 0, 0, 1, 0, 1, 1, 1]/2 - 0.25;
            wiggles = unique(wiggles.' .* wiggle_option, 'rows').';
            N_wiggles = size(wiggles, 2);
            wiggle_descriptors(N_wiggles) = no_wiggle; % pre-allocate memory
            for w_i=1:N_wiggles %calculate all wiggled coordinates
                wiggle_descriptors(w_i) = obj.wiggle_descriptor(wiggles(:, w_i));
            end
        end
        
        function wd = wiggle_descriptor(obj, dir)
            %
            % Construct shifted coordinates and phase ramps for a wiggle
            % steps. Pre-scale Fourier coordinates to optimize the
            % propagation functions a bit
            %
            wd = struct;
            % construct coordinates, shift quarter of a pixel when wiggling
            wd.pxe = obj.data_array((obj.grid.px_range - obj.grid.dpx * dir(2))/sqrt(obj.epsilon));
            wd.pye = obj.data_array((obj.grid.py_range - obj.grid.dpy * dir(1))/sqrt(obj.epsilon));
            wd.pze = obj.data_array((obj.grid.pz_range - obj.grid.dpz * dir(3))/sqrt(obj.epsilon));
            % construct phase gradients to compensate for the pixel
            % shift
            wd.gx = obj.data_array(exp(2.0i * pi * dir(2) * obj.grid.x_range / obj.grid.dx / length(obj.grid.x_range)));% - floor(obj.grid.N(2)/2));
            wd.gy = obj.data_array(exp(2.0i * pi * dir(1) * obj.grid.y_range / obj.grid.dx / length(obj.grid.y_range)));%
            wd.gz = obj.data_array(exp(2.0i * pi * dir(3) * obj.grid.z_range / obj.grid.dx / length(obj.grid.z_range)));%
        end
        
        function Ecrop = crop_field(obj,E)
            % Removes the boundary layer from the simulated field by
            % cropping field dataset
            Ecrop = E(obj.output_roi(1,1):obj.output_roi(2,1),...
                      obj.output_roi(1,2):obj.output_roi(2,2),...
                      obj.output_roi(1,3):obj.output_roi(2,3),...
                      obj.output_roi(1,4):obj.output_roi(2,4));
        end
        
        function rel_error = calculate_rel_error(obj, state)
            % Calculates relative error in final field by computing the
            % residual field: E - (GVE + GS).  dE_1 = GS, dE_k = GVdE_(k-1)
            % notes: 
            % no wiggle applied here so added errors are expected near boundaries. 
            % what to do with close zero-valued pixels?
            % todo: thoroughly test function 
            Etmp = -1.0i*obj.epsilon*obj.gamma.*state.E;         % E = VE
            Etmp = obj.propagate(Etmp,obj.no_wiggle)/obj.epsilon;% E = G'E/eps
            Etmp = state.source.add_to(Etmp,1.0i / obj.epsilon); % E = E + GS
                        
            rel_error = mean(abs(state.E(:) - Etmp(:)).^2)./mean(abs(state.E(:)).^2);
            rel_error = gather(rel_error);
        end
    end
end


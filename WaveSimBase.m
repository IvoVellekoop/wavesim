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
        gamma;     % potential array used in simulation (uses multiple arrays when medium_wiggle is enabled)
        epsilon;   % convergence parameter
        k02e;      % precomputed constants (pre-divided by sqrt epsilon)
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
        n_media;    
        % number of potential maps stored in memory (multiple submedia are 
        % used for anti-aliasing algorithm)
        
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
            
            % temporary fix for change in syntax (wiggle is now callled boundary_wiggle)
            if isfield(options,'wiggle')
                options.boundary_wiggle = options.wiggle;
            end
            % call simulation constructor
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
            obj.n_media = sample.n_media;
            obj.gamma = cell(obj.n_media,1);
            for i_medium = 1:sample.n_media                
                V = sample.e_r{i_medium}*k0^2-k0c^2 - 1.0i*obj.epsilon; % calculate potential map                   
                V = Medium.apply_filters(V,sample.filters); % apply sample filters to edge of potential map
                obj.gamma{i_medium} = obj.data_array(1.0i / obj.epsilon * V); % calculate gamma and convert to desired data type
            end

            %% convert to single or double precision, and put on gpu if
            % needed. 
            obj.k02e = obj.data_array(k0c^2/obj.epsilon + 1.0i);
            obj.epsilon = obj.data_array(obj.epsilon);      

            %% calculate wiggle descriptors
            [obj.wiggles, obj.boundary_wiggle, obj.medium_wiggle] = obj.compute_wiggles();    
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
            state.dE = obj.data_array([], obj.N);
            
            %% calculate number of wiggles
            Nwiggles = numel(obj.wiggles);             % total number of wiggles
            Nwiggles_per_medium= Nwiggles/obj.n_media; % number of boundary wiggles performed on a single submedium
            
            %% simulation iterations
            while state.has_next
                % select correct wiggle and medium number based on iteration number
                i_wiggle = mod(state.it-1, Nwiggles)+1; % wiggle counter
                i_medium = ceil(i_wiggle/Nwiggles_per_medium); % medium counter
                wig = obj.wiggles{i_wiggle};               
                
                % add source term (first Nwiggles iterations only) 
                if state.it <= Nwiggles % During the first few iterations: add source term
                    Etmp = state.source.add_to(state.dE, 1.0i / obj.epsilon / Nwiggles);
                else
                    Etmp = state.dE;
                end               
                
                % main computations steps
                Etmp = obj.propagate(Etmp, wig);
                state.dE = obj.mix(state.dE, Etmp, obj.gamma{i_medium}, wig);
                state.E = state.E + state.dE;

                % check if algorithm has to be terminated
                if state.calculate_energy
                    state.last_step_energy = Simulation.energy( obj.crop_field(state.dE) );                  
                end
                
                can_terminate = mod(state.it, Nwiggles) == 0; %only stop after multiple of Nwiggles iterations
                state = next(obj, state, can_terminate);
            end
            
            % Final steps
            % divide field by gamma to convert E' -> E
            wig = obj.wiggles{1};
            state.E = obj.wiggle_transform(state.E, wig.gpx, wig.gpy, wig.gpz);
            state.E = state.E ./ obj.gamma{1}; % this step will introduce aliasing artefacts
            state.E = obj.wiggle_transform(state.E, conj(wig.gpx), conj(wig.gpy), conj(wig.gpz));
            
            % crop field to match roi size
            state.E = obj.crop_field(state.E);
        end
    end
    methods(Access=protected)
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
        
        function dE = mix(obj, dE, Eprop, gamma, wig)
            % applies medium potential to previous difference field and the
            % propagated field and combines the two fields
            
            % transform the k-space wiggle phase ramp if medium_wiggle is
            % enabled.
            dE = obj.wiggle_transform(dE, wig.gpx, wig.gpy, wig.gpz);
            Eprop = obj.wiggle_transform(Eprop, wig.gpx, wig.gpy, wig.gpz);
            
            % mixes two (wiggled) fields
            if obj.gpu_enabled
                dE = arrayfun(@f_mix, dE, Eprop, gamma);
            else
                dE = f_mix(dE, Eprop, gamma);
            end
            
            % reverses k-space wiggle phase ramp of the combined field
            dE = obj.wiggle_transform(dE, conj(wig.gpx), conj(wig.gpy), conj(wig.gpz));
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
        function [wiggle_descriptors,boundary_wiggle,medium_wiggle] = compute_wiggles(obj) 
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

            if medium_wiggle == true                
                medium_wiggle = obj.N(1:3)' > 1; % medium wiggle disabled in directions where gridsize is 1
            elseif medium_wiggle == false
                medium_wiggle = false(3,1);
            elseif numel(medium_wiggle) < 3
                medium_wiggle(end+1:3) = false;
            end

            % determine wiggle directions
            wiggle_flags = [boundary_wiggle;medium_wiggle];
            wiggle_set = wiggle_perm(wiggle_flags);
            
            % calculate phase ramps and coordinates for every
            % boundary_wiggle and medium_wiggle
            Nwiggles = size(wiggle_set,2);
            wiggle_descriptors = cell(Nwiggles,1);  % pre-allocate memory
            for w_i=1:Nwiggles 
                wiggle_descriptors{w_i} = obj.wiggle_descriptor(wiggle_set(:,w_i));                
            end
        end
        
        function wd = wiggle_descriptor(obj, wig)
            % Constructs shifted coordinates and phase ramps for a wiggle
            % steps. (wig: logical vector 6x1 [y;x;z;ky;kx;kz])
            wd = struct;
            % construct coordinates, shift quarter of a pixel when wiggling
            % (required for calculating wiggled Green's function). Pre-scale 
            % Fourier coordinates to optimize the propagation functions a bit
            wd.pxe = obj.data_array((obj.grid.px_range - obj.grid.dpx * wig(2)/4)/sqrt(obj.epsilon));
            wd.pye = obj.data_array((obj.grid.py_range - obj.grid.dpy * wig(1)/4)/sqrt(obj.epsilon));
            wd.pze = obj.data_array((obj.grid.pz_range - obj.grid.dpz * wig(3)/4)/sqrt(obj.epsilon));
            
            % construct real space phase gradients to compensate for the
            % pixel shift in k_space (required for boundary_wiggle)
            wd.gx = obj.data_array(exp(1.0i * pi/2 * wig(2) * obj.grid.x_range / obj.grid.dx / length(obj.grid.x_range)));
            wd.gy = obj.data_array(exp(1.0i * pi/2 * wig(1) * obj.grid.y_range / obj.grid.dx / length(obj.grid.y_range)));
            wd.gz = obj.data_array(exp(1.0i * pi/2 * wig(3) * obj.grid.z_range / obj.grid.dx / length(obj.grid.z_range)));
            
            % construct k-space phase gradients to compensate for the
            % pixel shift in real space (required for medium_wiggle)
            wd.gpx = obj.data_array(exp(1.0i * pi/2 * wig(5) * obj.grid.px_range / obj.grid.dpx / length(obj.grid.px_range)));
            wd.gpy = obj.data_array(exp(1.0i * pi/2 * wig(4) * obj.grid.py_range / obj.grid.dpy / length(obj.grid.py_range)));
            wd.gpz = obj.data_array(exp(1.0i * pi/2 * wig(6) * obj.grid.pz_range / obj.grid.dpz / length(obj.grid.pz_range)));      
        end
        
        function E = fft_wiggle(obj, E, wig)
            % Modified Fourier transform: additionally applies wiggle phase
            % ramps in real-space. Todo: combine with wiggle_transform
            if obj.gpu_enabled
                E = arrayfun(@f_wiggle, E, wig.gx, wig.gy, wig.gz);
            else
                E = f_wiggle(E, wig.gx, wig.gy, wig.gz);              
            end
            E = fftn(E);
        end
        
        function E = ifft_wiggle(obj, E, wig)
            % Modified inverse Fourier transform: additionally reverses wiggle
            % phase ramps in real-space. Todo: combine with wiggle_transform
            E = ifftn(E);
            if obj.gpu_enabled              
                E = arrayfun(@f_wiggle, E, conj(wig.gx), conj(wig.gy), conj(wig.gz));
            else
                E = f_wiggle(E, conj(wig.gx), conj(wig.gy), conj(wig.gz));
            end
        end
        
        function E = wiggle_transform(obj, E, gpx, gpy, gpz)
            % Transforms k-space wiggle phase ramp of a field in real-space
            % Is required by the anti-aliasing algorithm. transform is only 
            % performed in the medium_wiggle direction that is enabled
            if obj.medium_wiggle(2) % px-direction
                E = ifft( fft(E,obj.grid.N(2),2) .* gpx , obj.grid.N(2), 2);
            end
            if obj.medium_wiggle(1) % py-direction
                E = ifft( fft(E,obj.grid.N(1),1) .* gpy , obj.grid.N(1), 1);
            end
            if obj.medium_wiggle(3) % pz-direction
                E = ifft( fft(E,obj.grid.N(3),3) .* gpz , obj.grid.N(3), 3);
            end
        end
    end
end


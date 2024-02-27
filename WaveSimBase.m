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
    % Ivo Vellekoop & Gerwin Osnabrugge 2016-2020
    
    properties
        k02e;      % precomputed constants (pre-divided by sqrt epsilon)
        gamma;     % potential array used in simulation
        epsilon;   % convergence parameter        
        filters;   % filters applied on the edge of the potential map       
        ACC = true;
        % 'true' indicates that the anti-cyclic convolution (or anti-wraparound) 
        % algorithm is  used on in all dimensions with non-zero boundary width. 
        % Note: zero-width boundaries are treated as periodic boundaries, 
        % and the anti-wraparound algorithm is disabled by default for these
        % boundaries. To override this default behavior, you can explicitly
        % pass a 3x1 logical vector (e.g. [true false false]) to indicate
        % which boundaries to 'wiggle'. 
        %% internal properties
        wiggles;
        % 'wiggling' is the multiplication of the field with a linear phase 
        % ramp this cell array containing all the phase ramps and the 
        % shifted Fourier space coordinates used for the anti-convolution 
        % convolution  
        epsilonmin = 3; 
        % minimum value to avoid divergence when simulating empty medium
        % (can perhaps be further optimized?)
    end
    methods(Abstract)
        propagate(obj); % function performing the propagation step
    end
    methods
        function obj = WaveSimBase(refractive_index, options)
            %% Constructs a wave simulation object
            % refractive_index   = refractive index map, may be complex, need not
            %                      be square. Can be 2-D or 3-D.
            % options.lambda = free space wavelength (same unit as pixel_size, e. g. um)
            % options.epsilon = convergence parameter (leave empty unless forcing a specific value)
            
            % temporary fix for change in syntax (wiggle is now callled ACC)
            if isfield(options,'wiggle')
                options.ACC = options.wiggle;
            end   
            
            % call simulation constructor
            obj@Simulation(refractive_index, options);
            fftw('planner','patient'); %optimize fft and ifft at first use
            
            %% determine wavenumbers
            % the optimal k_0 follows is given by n_center (see Medium)
            k0 = 2*pi/obj.lambda;                   % wavenumber in empty medium
            k0c = sqrt(obj.sample.e_r_center) * k0; % wavenumber for center value refractive index
            
            %% Determine epsilon and gamma
            % determine optimal epsilon (if value is not given in options)
            if ~isfield(options, 'epsilon') % setting epsilon in options forces a specific value, may not converge
                obj.epsilon = obj.calculate_epsilon(obj.sample.e_r,k0,k0c);
            end
            obj.iterations_per_cycle = obj.lambda /(2*k0c/obj.epsilon); %divide wavelength by pseudo-propagation length
            
            % calculate potential map and apply filters
            V = obj.sample.e_r*k0^2-k0c^2 - 1.0i*obj.epsilon; % calculate potential map
            [V,obj.filters] = obj.apply_edge_filters(V, options); % apply filters to edge of potential map
            
            % calculate potential map with background gain (gamma)
            obj.gamma = obj.data_array(1.0i / obj.epsilon * V);
            
            %% convert properties to single or double precision, and put on 
            % gpu if needed.                  
            obj.k02e = obj.data_array(k0c^2/obj.epsilon + 1.0i);
            obj.epsilon = obj.data_array(obj.epsilon);
            
            %% calculate wiggle descriptors
            [obj.wiggles, obj.ACC] = obj.compute_wiggles();
        end
        
        function state = run_algorithm(obj, state)
            % paper:
            %   E = [1+M+M^2+...] γ G S
            % iterate (with modification https://arxiv.org/abs/2207.14222)
            %   E_{0}   = 0
            %   E_{k+1} = E_k + α [(γ G V - γ) E_k + γ G S]
            %           = M E_k + α γ G S
            %   M = α (γ G V - γ) + I
            %
            % we now substitute V = -i ε γ and define G' =  ε G so that:
            %   M = -i α γ G' γ - γ + I

            % Now instead accumulate all in one buffer. dE_{k} are the terms in the Born series:
            %   dE_{k} = M^k α γ G S
            %
            % or, recursively:
            %   dE_{1} = α γ G S = i α γ / ε G' S
            %   dE_{k} = M dE_{k-1} = [-i α γ G' γ - γ + I] dE_{k-1}
            %
            % and we accumulate all in one buffer to find E
            %   E = dE_{1} + dE_{2} + ...
            %
            % A further optimization can be done by replacing dE = dE' / γ
            %
            %   dE'_{1} = i α γ² / ε G' S
            %   dE'_{k} = [-i α γ² G' - γ + I] dE'_{k-1}
            %           =  -i α γ² G' dE'_{k-1}  + (1-α γ) dE'_{k-1}
            %
            % Which simplifies to
            %
            %   dE'_{0} = 0
            %   dE'_{1} = -i α γ² G' [dE'_{k-1} + i/ε S]  + (1-α γ) dE'_{k-1}
            %   dE'_{k} = -i α γ² G'  dE'_{k-1}           + (1-α γ) dE'_{k-1}
            
            %   These iterations are implemented as:
            %   1. a propagation step:   Eprop   =  G' [E_diff + i/ε S] {only add S in 1st iteration}
            %   2. a mixing step:        E_diff  => (1-α γ) E_diff - i α γ² Eprop
            %   3. an accumulation step: E       => E + E_diff
            %
            %   And, after all iterations:
            %       E       => E / \gamma
            
            if(obj.gpu_enabled && obj.usemex )
                if(state.max_iterations > 0 && state.max_iterations < inf)
                   maxIter = state.max_iterations;
                else
                   maxIter = 1000;                 
                end
                
                callbackString = 'NoCallback';
                if(strcmp(func2str(obj.callback), 'Simulation.energy_added_disp_callback'))
                    callbackString = 'EnergyAddedDisp';
                elseif(strcmp(func2str(obj.callback), 'Simulation.energy_added_callback'))
                    callbackString = 'EnergyAdded';             
                end
                
                k02eSingle = single(gather(obj.k02e));
                epsilonSingle = single(gather(obj.epsilon));
                Nwiggles = numel(obj.wiggles);
                maxIter = double(gather(maxIter));

                [state.E, state.diff_energy] = RunWaveSimAlgorithmMex(obj.N, obj.gamma, Nwiggles, obj.wiggles, state.source, ...
                epsilonSingle, k02eSingle, obj.roi, maxIter, obj.energy_threshold, callbackString, obj.callback_interval);
                
                %% Final computional steps
                % divide field by gamma to convert E' -> E
                state.E = state.E ./ obj.gamma; 

                % crop field to match roi size
                state.E = obj.crop_field(state.E);
                state.it = length(state.diff_energy) +1;
                state.converged = state.it < maxIter;
            else
                %% Allocate memory for calculations
                state.E = obj.data_array([], obj.N);
                state.dE = obj.data_array([], obj.N);
                
                %% calculate number of different wiggles
                Nwiggles = numel(obj.wiggles);                    % total number of wiggles
    
                %% simulation iterations
                while state.has_next
                    % select correct wiggle and medium number based on iteration number
                    i_wiggle = mod(state.it-1, Nwiggles)+1; % wiggle counter
                    wig = obj.wiggles{i_wiggle};
                    
                    % add source term (first Nwiggles iterations only)
                    if state.it <= Nwiggles % During the first few iterations: add source term
                        Etmp = state.source.add_to(state.dE, 1.0i / obj.epsilon / Nwiggles);
                    else
                        Etmp = state.dE;
                    end
                    
                    % main computation steps
                    Etmp = obj.propagate(Etmp, wig);
                    state.dE = obj.mix(state.dE, Etmp, obj.gamma);
                    state.E = state.E + state.dE;
                    
                    % check if algorithm has to be terminated
                    if state.calculate_energy
                        state.last_step_energy = Simulation.energy( obj.crop_field(state.dE) );
                    end
                    
                    can_terminate = mod(state.it, Nwiggles) == 0; %only stop after multiple of Nwiggles iterations
                    state = next(obj, state, can_terminate);
                end
                
                %% Final computional steps
                % divide field by gamma to convert E' -> E
                state.E = state.E ./ obj.gamma; 
                
                % crop field to match roi size
                state.E = obj.crop_field(state.E);
            end
        end
    end
    methods(Access=protected)
        function epsilon = calculate_epsilon(obj, e_r, k0, k0c)
            % function used to calculate the optimal convergence parameter
            % epsilon used by wavesim based on maximum real value of
            % potential map
            
            % determine largest absolute potential fluctuation of all submedia
            % to guarantee convergence
            V_tot = e_r*k0^2-k0c^2;
            [Vabs_max, max_index] = max(abs(V_tot(:)));
            
            % For optimal convergence, reduce scaling with a factor of
            % 0.95 (see https://arxiv.org/abs/2207.14222)
            Vabs_max = Vabs_max / 0.95;
            epsilon = max(Vabs_max, obj.epsilonmin); % minimum value epsilonmin to avoid divergence when simulating empty medium
        end
        
        function [V,filters] = apply_edge_filters(obj, V, options)
            % filters the edges of the potential map V to reduce
            % reflections at the boundaries (is used for the ARL
            % boundaries)
            
            % construct edge filters
            filters = cell(3,1);
            
            % V is only filtered when ARL boundary type is used
            if (~strcmp(obj.sample.boundary_type, 'ARL') || all(obj.sample.Bl==0))
                return;
            end
            
            for dim=1:3
                bl = obj.sample.Bl(dim); %width of added boundary
                br = obj.sample.Br(dim);
                roi_size = obj.grid.N(dim) - bl - br;
                if bl > 0                
                    L = options.boundary_widths(dim); % width of boundary layer
                    if all(obj.grid.N > 1) % use linear window designed for 3D simulations
                        window = @(B) ((1:B)-0.21).'/(B+0.66);                      
                    else % use nuttall window for 1D and 2D simulations
                        win = @(nutall,B) nutall(1:B);
                        window = @(B) win(nuttallwin(2*L-1),B);
                    end
                    
                    % add a zero to the window if number of grid points is even 
                    if br == bl
                        smoothstep = [0; window(L-1)];
                    else % for odd number of grid points zeros are already added
                        smoothstep = window(L);
                    end
                    
                    % construct filter
                    filt = [zeros(bl-L, 1); smoothstep; ones(roi_size, 1); flipud(smoothstep); zeros(br-L, 1)];
                    filters{dim} = reshape(filt, circshift([1, 1, length(filt)], [0,dim]));
                    
                    % apply filter to potential map
                    V = V .* filters{dim};
                end
            end
        end
        
        function dE = mix(obj, dE, Eprop, gamma)
            % applies medium potential to previous difference field and the
            % propagated field and combines the two fields
            if obj.gpu_enabled
                dE = arrayfun(@f_mix, dE, Eprop, gamma);
            else
                dE = f_mix(dE, Eprop, gamma);
            end
        end
       
        
        %%% Wiggle methods (move to separate class?)
        function [wiggle_descriptors,ACC] = compute_wiggles(obj)
            % Decides which borders to wiggle and returns wiggled coordinates
            % for those borders
            
            % initialize boundary wiggle
            % true -> wiggle all non-periodic boundaries
            % [true, false, true] -> wiggle only 1st and 3rd dimension
            % false -> same as [false, false, false]
            ACC = obj.ACC(:);
            if numel(ACC) == 1
                if ACC == true
                    ACC = ~obj.grid.periodic(:);
                else
                    ACC = false(3,1);
                end
            end
            
            % determine all permutations of wiggle directions
            wiggle_set = wiggle_perm(ACC);
            
            % calculate phase ramps and coordinates for every different wiggle
            Nwiggles = size(wiggle_set,2);
            wiggle_descriptors = cell(Nwiggles,1);  % pre-allocate memory
            for w_i=1:Nwiggles
                wiggle_descriptors{w_i} = obj.wiggle_descriptor(wiggle_set(:,w_i));
            end
        end
        
        function wd = wiggle_descriptor(obj, wig)
            % Constructs shifted coordinates and phase ramps for a wiggle
            % steps. (wig: logical vector 3x1 [y;x;z])
            wd = struct;
            % construct coordinates, shift quarter of a pixel when wiggling
            % (required for calculating anti-cyclic Green's function). Pre-scale
            % Fourier coordinates to optimize the propagation functions a bit
            wd.pxe = obj.data_array((obj.grid.px_range - obj.grid.dpx * wig(2)/4)/sqrt(obj.epsilon));
            wd.pye = obj.data_array((obj.grid.py_range - obj.grid.dpy * wig(1)/4)/sqrt(obj.epsilon));
            wd.pze = obj.data_array((obj.grid.pz_range - obj.grid.dpz * wig(3)/4)/sqrt(obj.epsilon));
            
            % construct real space phase gradients to compensate for the
            % pixel shift in k_space (required for the ACC algorithm)
            wd.gx = obj.data_array(complex(exp(1.0i * pi/2 * wig(2) * obj.grid.x_range / obj.grid.dx / length(obj.grid.x_range))));
            wd.gy = obj.data_array(complex(exp(1.0i * pi/2 * wig(1) * obj.grid.y_range / obj.grid.dx / length(obj.grid.y_range))));
            wd.gz = obj.data_array(complex(exp(1.0i * pi/2 * wig(3) * obj.grid.z_range / obj.grid.dx / length(obj.grid.z_range))));
        end       
    end
end


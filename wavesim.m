classdef wavesim < simulation
    %Simulation of the 2-D or 3-D scalar wave equation using a Born series approach
    % Ivo M. Vellekoop 2014
    
    properties
        gamma;   % potential array used in simulation
        k;   % wave number for g0
        epsilon; %convergence parameter
        g0_k; % bare Green's function used in simulation, premultiplied by epsilon/i
        %% diagnostics and feedback
        epsilonmin; %minimum value of epsilon for which convergence is guaranteed (equal to epsilon, unless a different value for epsilon was forced)
        filters; %information for fft edge processing
    end
    
    methods
        function obj=wavesim(sample, options)
            %% Constructs a wave simulation object
            %	sample = SampleMedium object
            %   options.lambda = free space wavelength (same unit as pixel_size, e. g. um)
            %   options.epsilon = convergence parameter (leave empty unless forcing a specific value)
            obj@simulation(sample, options);
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use
            
            %% Determine k_0 to minimize epsilon
            % for now, we only vary the real part of k_0 and choose its
            % imaginary part as 0. Therefore, the optimal k_0
            % follows is given by n_center (see SampleMedium)
            %% Determine epsilon
            k00 = 2*pi/options.lambda;
            obj.k = sqrt(sample.e_r_center) * k00;
            
            % First construct V without epsilon
            V = sample.e_r*k00^2-obj.k^2;
            obj.epsilonmin = max(abs(V(:)));
            obj.epsilonmin = max(obj.epsilonmin, 3); %%minimum value to avoid divergence when simulating empty medium
            if isfield(options,'epsilon')
                obj.epsilon = options.epsilon*k00^2; %explicitly setting epsilon forces a specific value, may not converge
            else
                obj.epsilon = obj.epsilonmin; %guaranteed convergence
            end
            obj.iterations_per_cycle = obj.lambda /(2*obj.k/obj.epsilon); %divide wavelength by pseudo-propagation length
            
            %% Potential map (V==k^2-k_0^2-1i*epsilon)
            V = V - 1.0i*obj.epsilon;
            for d=1:3
                if ~isempty(sample.filters{d})
                    V = V .* sample.filters{d};
                end
            end
            obj.gamma = 1.0i / obj.epsilon * V;
            
            %% Calculate Green function for k_red (reduced k vector: k_red^2 = k_0^2 + 1.0i*epsilon)
            px = sample.grid.px_range;
            py = sample.grid.py_range;
            pz = sample.grid.pz_range;            
            obj.g0_k = -1.0i * obj.epsilon./(px.^2+py.^2+pz.^2-obj.k^2 - 1.0i*obj.epsilon);
            
            %convert to single our double precision, and put on gpu if
            %needed
            obj.gamma = obj.data_array(obj.gamma);
            obj.g0_k = obj.data_array(obj.g0_k);
            obj.epsilon = obj.data_array(obj.epsilon);
            obj.filters = sample.filters; %todo: remove this line
        end
        
        function state = run_algorithm(obj, state)
            % paper: 
            %   E = [1+M+M^2+...] \gamma G S
            % iterate:
            %   E_{0}   = 0
            %   E_{k+1} = M E_k + \gamma G S
            %   M = \gamma G V - \gamma + 1
            %   \gamma = i/epsilon V
            %
            % original implementation (equivalent):
            %   Ediff = i/epsilon V [E_k - G (V E_k + S)]
            %   E_{k+1} = E_{k} - Ediff
            %
            % optimized (no need to add source every step)
            % iterate:
            %   Ediff_{0}   = \gamma G S
            %   E_{0}       = \gamma G S
            %
            %   Ediff_{k+1} = M Ediff_{k} = (\gamma G V - \gamma + 1) E_diff_{k}
            %   E_{k+1}     = E_{k} + Ediff_{k+1}
            %   
            %   implementation:
            %   Ediff_{k+1} = M Ediff_{k} = (1 - \gamma) E_diff_{k} + \gamma G' \gamma E_diff_{k}
            %   with G' = G epsilon / i so that G'\gamma = G V
            
            
            %% Allocate memory for calculations
            % start iteration with the \gamma G S term:
            Ediff = simulation.add_sources(state, data_array(obj), obj.roi);
            Ediff = obj.gamma .* ifftn(1.0i / obj.epsilon * obj.g0_k .* fftn(Ediff)); %gamma G S
            state.E = Ediff;
            
            %% simulation iterations
            while state.has_next
                Ediff = (1-obj.gamma) .* Ediff + obj.gamma .* ifftn(obj.g0_k .* fftn(obj.gamma .* Ediff));
           
                if state.calculate_energy
                   state.last_step_energy = simulation.energy(Ediff, obj.roi);
                end
                
                state.E = state.E + Ediff;
                state = next(obj, state);
            end
        end
    end
end
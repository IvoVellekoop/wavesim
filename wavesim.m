classdef wavesim < simulation
    %Simulation of the 2-D or 3-D scalar wave equation using a Born series approach
    % Ivo M. Vellekoop 2014
    
    properties
        V;   % potential array used in simulation
        Rp;   % prefactor for resummed approach
        k;   % wave number for g0
        epsilon; %convergence parameter
        g0_k; % bare Green's function used in simulation (center subtracted in resummed approach)
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
            obj.V = sample.e_r*k00^2-obj.k^2;
            obj.epsilonmin = max(abs(obj.V(:)));
            obj.epsilonmin = max(obj.epsilonmin, 3); %%minimum value to avoid divergence when simulating empty medium
            if isfield(options,'epsilon')
                obj.epsilon = options.epsilon*k00^2; %explicitly setting epsilon forces a specific value, may not converge
            else
                obj.epsilon = obj.epsilonmin; %guaranteed convergence
            end
            obj.iterations_per_cycle = obj.lambda /(2*obj.k/obj.epsilon); %divide wavelength by pseudo-propagation length
            
            %% Potential map (V==k^2-k_0^2-1i*epsilon)
            obj.V = obj.V - 1.0i*obj.epsilon;
            for d=1:3
                if ~isempty(sample.filters{d})
                    obj.V = obj.V .* sample.filters{d};
                end
            end
            
            %% Calculate Green function for k_red (reduced k vector: k_red^2 = k_0^2 + 1.0i*epsilon)
            obj.g0_k = 1./(p2(sample.grid) - obj.k^2 - 1.0i*obj.epsilon);
            
            %%%% resummation: move center to origin and scale to unit circle to get \Phi
            obj.g0_k = - 2.0i * obj.epsilon * obj.g0_k - 1; %\Phi
            obj.Rp = (1.0i * obj.V) ./ (2 * obj.epsilon - 1.0i * obj.V ); 
            
            %convert to single or double precision, and put on gpu if
            %needed
            obj.V = obj.data_array(obj.V);
            obj.g0_k = obj.data_array(obj.g0_k);
            obj.epsilon = obj.data_array(obj.epsilon);
            obj.filters = sample.filters; %todo: remove this line
        end
        
        function state = run_algorithm(obj, state)
            %% Allocate memory for calculations
            Esrc = simulation.add_at(data_array(obj), state.source, state.source_pos);
            G = 1.0i / (2 * obj.epsilon) * (1 + obj.g0_k);
            Esrc = - obj.V .* obj.V / obj.epsilon^2 .* ifftn(G .* fftn(Esrc));
            state.E = data_array(obj);    
            
            %% simulation iterations
            while state.has_next
           %original:
           %     Etmp = simulation.add_at(obj.V.* state.E, state.source, state.source_pos);
           %     Ediff = -(1.0i/obj.epsilon*obj.V) .* (state.E-ifftn(obj.g0_k .* fftn(Etmp)));
           
           %resummed:
                state.E = obj.Rp .* ifftn(obj.g0_k .* fftn(state.E + Esrc));
           
                if state.calculate_energy
                   state.last_step_energy = 100;%simulation.energy(Ediff(obj.roi{1}, obj.roi{2}, obj.roi{3}));
                end
                
                %state.E = state.E + Ediff;
                state = next(obj, state);
            end
            state.E = (state.E + obj.Rp .* Esrc) ./ (obj.V * 1.0i / obj.epsilon) * 2;
        end
    end
end
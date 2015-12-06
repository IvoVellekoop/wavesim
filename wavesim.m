classdef wavesim < simulation
    %Simulation of the 2-D wave equation using a Born series approach
    % Ivo M. Vellekoop 2014
    
    properties
        V;   % potential array used in simulation
        differential_mode = false; %when set to 'true', only the differential field for each iteration is calculated: the fields are not added to get a solution to the wave equation (used for debugging)
        energy_threshold = 1E-20; %the simulation is terminated when the added energy between two iterations is lower than 'energy_threshold'. Default 1E-9
        %%internal
        k; % wave number
        epsilon;
        g0_k; % bare Green's function used in simulation
        %% diagnostics and feedback
        epsilonmin; %minimum value of epsilon for which convergence is guaranteed (equal to epsilon, unless a different value for epsilon was forced)
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
            obj.epsilonmin = max(obj.epsilonmin, 1E-3); %%minimum value to avoid divergence when simulating empty medium
            if isfield(options,'epsilon')
                obj.epsilon = options.epsilon*k00^2; %force a specific value, may not converge
            else
                obj.epsilon = obj.epsilonmin; %guaranteed convergence
            end
            
            %% Potential map (V==k^2-k_0^2-1i*epsilon)
            obj.V = obj.V - 1.0i*obj.epsilon;
            
            %% Calculate Green function for k_red (reduced k vector: k_red^2 = k_0^2 + 1.0i*epsilon)
            f_g0_k = @(px, py) 1./(px.^2+py.^2-(obj.k^2 + 1.0i*obj.epsilon));
            obj.g0_k = bsxfun(f_g0_k, sample.grid.px_range, sample.grid.py_range);
            
            if obj.gpuEnabled
                obj.V = gpuArray(obj.V);
                obj.g0_k = gpuArray(obj.V);
            end
        end
        
        function state = run_algorithm(obj, state)
            %% Allocate memory for calculations
            state.E = data_array(obj);
            
            %% simulation iterations
            while state.has_next
                Ediff = (1.0i/obj.epsilon*obj.V) .* (state.E-ifft2(obj.g0_k .* fft2(obj.V.*state.E+state.source)));
                if state.calculate_energy
                   state.last_step_energy = simulation.energy(Ediff(obj.roi{1}, obj.roi{2}));
                end
                state.E = state.E - Ediff;
                state = next(obj, state);
            end
        end
    end
%         function options = readout_input(options)
%             % function reading out all given constructor options and filling
%             % in default values for missing properties
%             if nargin == 0
%                 options = struct;
%             end
%             if ~isfield(options,'lambda')
%                 options.lambda = 1; % in um
%             end
%             % stopping criterion
%             if ~isfield(options,'energy_threshold')
%                 options.energy_threshold = 1E-20;
%             end
%             % callback function is called every N frames
%             if ~isfield(options,'callback_interval')
%                 options.callback_interval = 500;
%             end
%         end
end
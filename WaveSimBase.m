classdef(Abstract) WaveSimBase < Simulation
    % Base class for the simulation of 2-D or 3-D wave equations using
    % the modified Born series approach
    % This class is overridden by wavesim (for scalar simulations)
    % and wavesimv (for vector simulations).
    %
    % The actual simulation is performed by calling the functions
    % 'mix' and 'propagate'. Only the 'propagate' function is different
    % for vector or scalar waves.
    %
    % Ivo M. Vellekoop 2014-2018
    
    properties
        gamma;   % potential array used in simulation
        k;       % wave number for g0
        epsilon; % convergence parameter
        mix;     % function handle to function performing the mixing step
        wiggle;  % 'true' indicates that the anti-wraparound algorithm is used
        %% diagnostics and feedback
        epsilonmin; %minimum value of epsilon for which convergence is guaranteed (equal to epsilon, unless a different value for epsilon was forced)
    
        % precomputed vectors and constants (pre-divided by sqrt epsilon)
        pxe;
        pye;
        pze;
        k02e;
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
            if isfield(options, 'epsilon')
                obj.epsilon = options.epsilon; %explicitly setting epsilon forces a specific value, may not converge
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
            
            %% calculate wiggle coefficients
            if obj.wiggle
                %wiggle left/right in the dimensions that are not periodic
                %for simplicity, we always wiggle in all three dimensions,
                %but we reduce the wiggle amplitude to 0 in the periodic
                %directions
                wiggles = ...
                    [0, 1, 0, 1, 0, 1, 0, 1;...
                     0, 0, 1, 1, 0, 0, 1, 1;...
                     0, 0, 0, 0, 1, 1, 1, 1]/2 - 0.25;
                wiggles(obj.grid.periodic, :) = 0;
            else
                wiggles = [0;0;0];
            end
            %
            % Construct shifted coordinates and phase ramps for all wiggle
            % steps. Pre-scale Fourier coordinates to optimize the
            % propagation functions a bit
            %
            for w_i=1:size(wiggles, 2)
                w = wiggles(:, w_i);
                c = struct;
                % construct coordinates, shift half a pixel when wiggling
                c.pxe = obj.data_array((obj.grid.px_range - obj.grid.dpx * w(2))/sqrt(obj.epsilon));
                c.pye = obj.data_array((obj.grid.py_range - obj.grid.dpy * w(1))/sqrt(obj.epsilon));
                c.pze = obj.data_array((obj.grid.pz_range - obj.grid.dpz * w(3))/sqrt(obj.epsilon));
                % construct phase gradients to compensate for the pixel
                % shift
                c.gx = obj.data_array(exp(2.0i * pi * w(2) * obj.grid.x_range / obj.grid.dx / length(obj.grid.x_range)));% - floor(obj.grid.N(2)/2));
                c.gy = obj.data_array(exp(2.0i * pi * w(1) * obj.grid.y_range / obj.grid.dx / length(obj.grid.y_range)));%
                c.gz = obj.data_array(exp(2.0i * pi * w(3) * obj.grid.z_range / obj.grid.dx / length(obj.grid.z_range)));%
                if w_i == 1
                    obj.wiggle = c;
                else
                    obj.wiggle(w_i) = c;
                end
            end
            
            %% convert to single or double precision, and put on gpu if
            % needed. 
            obj.k02e = obj.data_array(obj.k^2/obj.epsilon + 1.0i);
            obj.gamma = obj.data_array(obj.gamma);
            obj.epsilon = obj.data_array(obj.epsilon);

            %%define function for the mixing step:
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
            %   Ediff_{0}   = S
            %   E_{0}       = 0
            %
            %   Ediff_{k+1} = M Ediff_{k} = (\gamma G V - \gamma + 1) E_diff_{k}
            %   E_{k+1}     = E_{k} + Ediff_{k+1}
            %   
            %   substitute V = -1.0i \epsilon \gamma
            %   Ediff_{k+1} = M Ediff_{k} = (1 - \gamma) E_diff_{k} - 1.0i \epsilon \gamma G \gamma E_diff_{k}
            %
            %   substitute G' = G epsilon  and Ediff' = \gamma Ediff
            %   (in the end we have to divide E' by \gamma to get E)
            %   Ediff'_{0}   = 0
            %   Ediff'_{k+1} = (1 - \gamma) E_diff'_{k} - 1.0i \gamma^2 G' E'_diff_{k}
            %
            %   These iterations are implemented as:
            %   1. a propagation step:   Eprop   = (G' E'_diff + 1.0i/epsilon S {only in 1st iteration})
            %   2. a mixing step:        E_diff' => (1-\gamma) E_diff' - 1.0i \gamma^2 Eprop  
            %   3. an accumulation step: E'      => E' + E_diff'
            %
            %   With finally a correction step: E    = E' / \gamma
            
            %% Allocate memory for calculations
            % start iteration with the source term (pre-multiplied by 1.0i/obj.epsilon. epsilon to compensate for pre-multiplication of G):
            % pre-multiply by gamma/epsilon
            state.E = 0;
            state.Ediff = obj.data_array([], obj.N);
            
            %% simulation iterations
            Nwiggle = size(obj.wiggle, 2);
            while state.has_next
                wigg = obj.wiggle(mod(state.it, Nwiggle) + 1); 
                if state.it <= size(obj.wiggle, 2)
                    state.Ediff = obj.mix(state.Ediff, obj.propagate(...
                        state.source.add_to(state.Ediff, 1.0i / obj.epsilon / Nwiggle)...
                        , wigg), obj.gamma);
                else
                    state.Ediff = obj.mix(state.Ediff, obj.propagate(state.Ediff, wigg), obj.gamma);
                end
                
                state.E = state.E + state.Ediff(...
                    obj.output_roi(1,1):obj.output_roi(2,1),...
                    obj.output_roi(1,2):obj.output_roi(2,2),...
                    obj.output_roi(1,3):obj.output_roi(2,3),...
                    obj.output_roi(1,4):obj.output_roi(2,4));
                
                if state.calculate_energy
                    state.last_step_energy = Simulation.energy(state.Ediff);
                end
                
                state = next(obj, state);
                if mod(state.it, Nwiggle) ~= 0
                    state.has_next = true; %only stop after multiple of 8 iterations
                end
            end
            state.E = state.E ./ obj.gamma(...
                obj.output_roi(1,1):obj.output_roi(2,1),...
                    obj.output_roi(1,2):obj.output_roi(2,2),...
                    obj.output_roi(1,3):obj.output_roi(2,3),...
                    obj.output_roi(1,4):obj.output_roi(2,4));
        end
    end
end


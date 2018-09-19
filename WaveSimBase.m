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
        wiggle = true; 
        no_wiggle; %holds coordinates for non-wiggled case
        % 'true' indicates that the anti-wraparound algorithm is used on
        % all edges with non-zero boundary width. Note: zero-width
        % boundaries are treated as periodic boundaries, and the
        % anti-wraparound algorithm is disabled by default for these
        % boundaries. To override this default behavior, you can explicitly
        % pass a 1x3 logical vector (e.g. [true false false]) to indicate
        % which boundaries to 'wiggle'
        %% diagnostics and feedback
        epsilonmin; %minimum value of epsilon for which convergence is guaranteed (equal to epsilon, unless a different value for epsilon was forced)
    
        % precomputed vectors and constants (pre-divided by sqrt epsilon)
        %pxe;
        %pye;
        %pze;
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
            % todo: move to Medium object
            V = V - 1.0i*obj.epsilon;
            for d=1:3
                if ~isempty(sample.filters{d})
                    V = V .* sample.filters{d};
                end
            end
            obj.gamma = 1.0i / obj.epsilon * V;
            
            %% calculate wiggle coefficients
            [obj.wiggle, obj.no_wiggle] = obj.compute_wiggles(obj.wiggle);
            
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
            % Now instead accumulate all in one buffer. dE_{k} are the terms in the Born series:
            %   dE_{1} = \gamma G S
            %   dE_{k} = M dE_{k-1}
            %   E = dE_{1} + dE_{2} + ...
            %
            %   substitute V = -1.0i \epsilon \gamma
            %   and G = G' / (-1.0i \epsilon)
            %
            %   dE_{1} = \gamma G' i/\epsilon S
            %   dE_{k} = (\gamma G' \gamma - \gamma + 1) dE_{k-1}
            %
            % A small optimization can be done by replacing dE = dE' / \gamma
            %
            %   dE'_{1} = \gamma^2 G' i/\epsilon S
            %   dE'_{k} = (\gamma^2 G' - \gamma + 1) dE'_{k-1}
            %           = \gamma^2 G' dE'_{k-1}  + (1-\gamma) dE'_{k-1}
            %
            % Which simplifies to
            %
            %   dE'_{0} = 0
            %   dE'_{1} = \gamma^2 G' [dE'_{k-1} + i/\epsilon S]  + (1-\gamma) dE'_{k-1}
            %   dE'_{k} = \gamma^2 G'  dE'_{k-1}                  + (1-\gamma) dE'_{k-1}
            
            %   These iterations are implemented as:
            %   1. a propagation step:   Eprop   = G' \gamma (E_diff + 1.0i/epsilon S) {only add S in 1st iteration}
            %   2. a mixing step:        E_diff  => (1-\gamma) E_diff - 1.0i \gamma^2 Eprop  
            %   3. an accumulation step: E       => E + E_diff
            %   
            %   And, after all iterations:
            %       E       => E / \gamma
            
            %% Allocate memory for calculations
            state.E = 0;
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
                
                % todo: test if it is faster to accumulate the full Ediff
                % and take the roi later.
                state.E = state.E + state.Ediff(...
                    obj.output_roi(1,1):obj.output_roi(2,1),...
                    obj.output_roi(1,2):obj.output_roi(2,2),...
                    obj.output_roi(1,3):obj.output_roi(2,3),...
                    obj.output_roi(1,4):obj.output_roi(2,4));
                
                if state.calculate_energy
                    state.last_step_energy = Simulation.energy(state.Ediff);
                end
                
                can_terminate = mod(state.it, Nwiggle) == 0; %only stop after multiple of Nwiggle iterations
                state = next(obj, state, can_terminate);                
            end
            
            state.E = state.E ./ obj.gamma(...
                    obj.output_roi(1,1):obj.output_roi(2,1),...
                    obj.output_roi(1,2):obj.output_roi(2,2),...
                    obj.output_roi(1,3):obj.output_roi(2,3),...
                    obj.output_roi(1,4):obj.output_roi(2,4));
        end
    end
    methods(Access=private)
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
                [0, 1, 0, 1, 0, 1, 0, 1;...
                 0, 0, 1, 1, 0, 0, 1, 1;...
                 0, 0, 0, 0, 1, 1, 1, 1]/2 - 0.25;
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
            % construct coordinates, shift half a pixel when wiggling
            wd.pxe = obj.data_array((obj.grid.px_range - obj.grid.dpx * dir(2))/sqrt(obj.epsilon));
            wd.pye = obj.data_array((obj.grid.py_range - obj.grid.dpy * dir(1))/sqrt(obj.epsilon));
            wd.pze = obj.data_array((obj.grid.pz_range - obj.grid.dpz * dir(3))/sqrt(obj.epsilon));
            % construct phase gradients to compensate for the pixel
            % shift
            wd.gx = obj.data_array(exp(2.0i * pi * dir(2) * obj.grid.x_range / obj.grid.dx / length(obj.grid.x_range)));% - floor(obj.grid.N(2)/2));
            wd.gy = obj.data_array(exp(2.0i * pi * dir(1) * obj.grid.y_range / obj.grid.dx / length(obj.grid.y_range)));%
            wd.gz = obj.data_array(exp(2.0i * pi * dir(3) * obj.grid.z_range / obj.grid.dx / length(obj.grid.z_range)));%
        end
    end
end


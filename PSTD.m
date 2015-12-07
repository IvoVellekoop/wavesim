classdef PSTD < simulation
    %Simulation of the 2-D wave equation using PSTD
    % Saroch Leedumrongwatthanakun 2015
    
    properties
        %optiond:
        dt_relative = 0.95; %Relative size of a time step compared to the convergence criterion for PSTD. Decrease for more accurate results (default=0.95)
        source_amplitude = @PSTD.default_source; %amplitude of the source as a function of time. Can be replaced by a different function by the user
        
        %internal
        c1;   % coefficients
        c2;
        c3;
        dt; % time step
        koperator; % k domain laplacian
        dtmax; %maximum time step at sutible pixel_size
    end
    
    methods
        function obj = PSTD(sample, options)
            %% Constructs a pseudo spectral time domain simulation object
            %	sample = SampleMedium object
            %   options.wavelength = free space wavelength (same unit as pixel_size, e. g. um)
            %   options.dt         = time step (leave empty unless forcing a specific value)
            %   options.time_duration = time duration of whole simulation
            obj@simulation(sample, options);
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use
            
            %% Calculate time step dt
            % The light speed is defined as 1 distance_unit / time_unit.
            %
            obj.dtmax = 2/sqrt(2)/pi*sample.grid.dx*sqrt(sample.e_r_min); %Stability condition (ref needed)
            obj.dt = obj.dt_relative * obj.dtmax;
            
            %% Initialize coefficients (could be optimized);
            sdt = imag(sample.e_r) * 2*pi/options.lambda * obj.dt;
            obj.c1 = (sdt-2)./(sdt+2);
            obj.c2 = 4./(sdt+2);
            obj.c3 = 2*obj.dt^2./real(sample.e_r)./(sdt+2);
            
            %% Calculate PSTD laplace operator
            f_laplace = @(px, py) -(px.^2+py.^2); %-k^2
            obj.koperator = (bsxfun(f_laplace, sample.grid.px_range, sample.grid.py_range));
        end
        
        function state = run_algorithm(obj, state)
            %% Allocate memory for calculations
            state.E = data_array(obj);
            E_prev = data_array(obj);
            A = 1; %source amplitude
            A_prev = 1; %previous source amplitude (used to calculate energy difference)
            %todo: gpuarray for c1,c2,c3 and koperator
            
            %% iterate algorithm
            while state.has_next
                %% Calculate source amplitude
                % First calculate time in number of optical cycles:
                % The light speed is defined as 1 distance_unit / time_unit.
                % Therefore, to convert iteration number to optical cycle:
                %  Ncycle = Nit * dt[time_unit] / lambda[distance_unit] 
                %  * c=1[distance_unit/time_unit]
                %
                Ncycle = state.it * obj.dt / obj.lambda;
                A = obj.source_amplitude(Ncycle);

                %% update fields
                %based on equation:
                % nabla^2 E - e_r d^2 E/dt^2 = -source
                %
                %approximate time derivative by centered finite difference
                %
                % nabla^2 E - e_r (E_next - 2*E + E_prev) / dt^2 = -source
                %
                % isolate E_next:
                % E_next = (nabla^2 E + source) * dt^2/e_r + 2*E - E_prev
                E_next = obj.c2.*state.E + obj.c1.*E_prev + obj.c3 .* (ifft2(obj.koperator.*fft2(state.E)) + A*state.source);
                
                if state.calculate_energy
                    phase_shift = exp(1.0i*(angle(A)-angle(A_prev))); %expected phase shift for single step (only works for CW source!!)
                    state.last_step_energy = wavesim.energy(E_next(obj.roi{1}, obj.roi{2}) - state.E(obj.roi{1}, obj.roi{2})*phase_shift);
                    if (A<0.5) %workaround: don't terminate when source is still spinning ups
                        state.last_step_energy = max(state.last_step_energy, state.threshold*2);
                    end;
                end
                %% update fields
                E_prev = state.E;
                A_prev = A;
                state.E = E_next;
                state = next(obj, state);
            end
            
            %finally, compensate for phase of source
            state.E = state.E / A;
        end
    end
    methods(Static)
        function A = default_source(Ncycle)
            %represents a CW source with a frequency corresponding to 3e8/lambda
            %the source is switched on slowly with a cosine profile (currently takes 100 cycles)
            %% calculate source amplitude
            A = exp(-2.0i*pi * Ncycle);
            %implement smooth start
            A = A .* (0.5-0.5*cos(pi*min(Ncycle/100.0,1.0)));
        end
    end
end




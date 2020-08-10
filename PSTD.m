classdef PSTD < Simulation
    % Simulation of the 2-D and 3-D wave equation using PsuedoSpectral Time
    % Domain
    % Saroch Leedumrongwatthanakun 2015
    
    properties
        %option:
        dt_relative = 0.95; %Relative size of a time step compared to the convergence criterion for PSTD. Decrease for more accurate results (default=0.95)
        source_amplitude = @PSTD.default_source; %amplitude of the source as a function of time. Can be replaced by a different function by the user
        
        %internal coefficients
        c1;
        c2;
        c3;
        dt; % time step
        omega; %2pi/lambda (c===1)
        koperator; % k domain laplacian
        dtmax; %maximum time step at sutible pixel_size
    end
    
    methods
        function obj = PSTD(refractive_index, options)
            warning('The PSTD code is only provided for comparison with wavesim, it is not optimized for real use.');
            %% Constructs a pseudo spectral time domain simulation object
            %  refractive_index   = refractive index map, may be complex, need not
            %                      be square. Can be 2-D or 3-D.
            %  options.wavelength = free space wavelength (same unit as pixel_size, e. g. um)
            %  options.dt         = time step (leave empty unless forcing a specific value)
            %  options.time_duration = time duration of whole simulation     
            obj@Simulation(refractive_index, options);
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use
            
            % only use first submedium (medium_wiggle not implemented for
            % PSTD)
            obj.sample.e_r = obj.sample.e_r{1};
            
            %% Calculate time step dt
            % The light speed is defined as 1 distance_unit / time_unit.
            %
            % determine dimensionality of the simulation for calculating
            % the maximum time step
            if min(obj.grid.N) == 1
                dimensions = 2;
            else
                dimensions = 3;
            end
            
            obj.dtmax = 2/sqrt(dimensions)/pi*obj.grid.dx*sqrt(obj.sample.e_r_min); %Stability condition
            obj.dt = obj.dt_relative * obj.dtmax;
            obj.iterations_per_cycle = obj.lambda / obj.dt; %lambda[distance] / dt[time] / c[distance/time]
            
            %% Initialize coefficients (could be optimized);
            obj.omega = 2*pi/obj.lambda; %wave speed c_0 = 1 distance unit / time unit by definition, so omega=k00
            c2dt = obj.dt^2./real(obj.sample.e_r); %(relative wave speed * dt) ^2
            sdt  = obj.dt*obj.omega*imag(obj.sample.e_r)./real(obj.sample.e_r); %sigma dt
            
            obj.c1 = (sdt-2)./(sdt+2);
            obj.c2 = 4./(sdt+2);
            obj.c3 = 2*c2dt./(sdt+2);
            
            %% Calculate PSTD laplace operator
            px = obj.grid.px_range;
            py = obj.grid.py_range;
            pz = obj.grid.pz_range;            
            obj.koperator = -px.^2-py.^2-pz.^2;
            obj.max_cycles = obj.max_cycles+100; %slow starting source
            
            if obj.gpu_enabled
                obj.c1 = gpuArray(obj.c1);
                obj.c2 = gpuArray(obj.c2);
                obj.c3 = gpuArray(obj.c3);
                obj.koperator = gpuArray(obj.koperator);
            end
        end
        
        function state = run_algorithm(obj, state)
            %% Allocate memory for calculations
            state.E = obj.data_array([], obj.N);
            E_prev = obj.data_array([], obj.N);
            A = 1; %source amplitude

            %% iterate algorithm
            while state.has_next
                %% Calculate source amplitude
                % First calculate time in number of optical cycles:
                % The light speed is defined as 1 distance_unit / time_unit.
                % Therefore, to convert iteration number to optical cycle:
                %  Ncycle = Nit * dt[time_unit] / lambda[distance_unit] 
                %  * c=1[distance_unit/time_unit]
                %
                Ncycle = state.it / obj.iterations_per_cycle;
                A_prev = A;
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
                Etmp = ifftn(obj.koperator.*fftn(state.E));
                Etmp = state.source.add_to(Etmp, A); 
                %simulation.add_at(ifftn(obj.koperator.*fftn(state.E)), A * state.source, state.source_pos);
                E_next = obj.c2.*state.E + obj.c1.*E_prev + obj.c3 .* Etmp;
                
                if state.calculate_energy
                    phase_shift = exp(1.0i*(angle(A)-angle(A_prev))); %expected phase shift for single step (only works for CW source!!)
                    state.last_step_energy = Simulation.energy(E_next - state.E * phase_shift, obj.roi) / abs(A)^2;
                end
                %% update fields
                E_prev = state.E;
                state.E = E_next;
                can_terminate = A > 0.9; %don't terminate when source is still spinning up
                state = next(obj, state, can_terminate);
            end
            
            %finally, crop to output_roi and compensate for phase of source
            A_next = obj.source_amplitude(state.it / obj.iterations_per_cycle);
            state.E = obj.crop_field(state.E) * exp(-1.0i*angle(A_next));
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




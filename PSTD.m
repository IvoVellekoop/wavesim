classdef PSTD < simulation
    %Simulation of the 2-D wave equation using PSTD
    % Saroch Leedumrongwatthanakun 2015
    
    properties
        c1;   % coefficients
        c2;
        c3;
        info; % diagnostics information on cover of field on medium
        dt; % time step
        number_of_time_steps; % a number of time step
        cover_time; % time of wave propagation covering roi
        koperator; % k domain laplacian
        dtmax; %maximum time step at sutible pixel_size
    end
    
    methods
        function obj = PSTD(sample, options)
            %% Constructs a pseudo spectral time domain simulation object
            % uses algorithm from 'Numerical Absorbing Boundary Conditions
            % Based on a Damped Wave Equation for Pseudospectral
            % Time-Domain Acoustic Simulations' Spa, Reche-Lopez, Hernandez
            % The Scientific World Journal 2014, 285945 (2014)
            %	sample = SampleMedium object
            %   options.wavelength = free space wavelength (same unit as pixel_size, e. g. um)
            %   options.dt         = time step (leave empty unless forcing a specific value)
            %   options.time_duration = time duration of whole simulation
            obj@simulation(sample, options);
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use
                        
            %% Stability condition (c=1)
            obj.dtmax = 2/sqrt(2)/pi*sample.grid.dx*sqrt(sample.e_r_min); %is this correct?
            if isfield(options,'dt')
                obj.dt = options.dt; %force a specific value, may not converge
            else
                obj.dt = 0.5*obj.dtmax; %guaranteed convergence
            end
            
            %% Initialize coefficients (could be optimized);
            sdt = imag(sample.e_r) * 2*pi/options.lambda * obj.dt;
            obj.c1 = (sdt-2)./(sdt+2);
            obj.c2 = 4./(sdt+2);
            obj.c3 = 2*obj.dt^2./real(sample.e_r)./(sdt+2); 
            
            %%  Calculate PSTD laplace operator
            f_laplace = @(px, py) -(px.^2+py.^2); %-k^2
            obj.koperator = (bsxfun(f_laplace, sample.grid.px_range, sample.grid.py_range));
        end
        
        function state = run_algorithm(obj, state)
            %% Allocate memory for calculations
            state.E = data_array(obj);
            E_prev = data_array(obj);
            %todo: gpuarray for c1,c2,c3 and koperator
            
            omega  = 2*pi/obj.lambda; %c==1
            
            while state.has_next
                %% calculate source amplitude
                osc = state.it * obj.dt * omega; 
                source_amplitude = exp(-1.0i * osc);
                if (osc < 100*pi)
                    source_amplitude = source_amplitude * (1-0.5*cos(osc/100));
                end
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
                E_next = obj.c2.*state.E + obj.c1.*E_prev + obj.c3 .* (ifft2(obj.koperator.*fft2(state.E)) + source_amplitude*state.source);
                
                if state.calculate_energy
                    phase_shift = exp(-1.0i*obj.dt*omega); %expected phase shift for single step
                    state.last_step_energy = wavesim.energy(E_next(obj.roi{1}, obj.roi{2}) - state.E(obj.roi{1}, obj.roi{2})*phase_shift);
                end
                %% update fields
                E_prev = state.E;
                state.E = E_next;
                state = next(obj, state);
            end
            
            %finally, compensate for phase of source
            state.E = state.E / source_amplitude;
        end 
    end
end
    
    


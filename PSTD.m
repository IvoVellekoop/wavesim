classdef PSTD
    %Simulation of the 2-D wave equation using PSTD
    % Saroch Leedumrongwatthanakun 2015
    
    properties
        c1;   % coefficients
        c2;
        c3;
        grid; %simgrid object
        roi; %position of simulation area with respect to padded array
        x_range;
        y_range;
        lambda;
        
        info; % diagnostics information on cover of field on medium
		gpuEnabled = false; % logical to check if simulation are ran on the GPU (default: false)
        callback = @PSTD.default_callback; %callback function that is called for showing the progress of the simulation. Default shows image of the real value of the field.
        callback_interval = 100; %the callback is called every 'callback_interval' steps. Default = 5
        dt; % time step
        it; %iteration
        time; % time consumption
        energy_threshold = 1E-9; %the simulation is terminated when the added energy between two iterations is lower than 'energy_threshold'. Default 1E-9
        max_iterations = 4E3; %or when 'max_iterations' is reached. Default 10000
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
            
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use
            options = PSTD.readout_input(options); %fill in default options
                        
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
            kOper2 = @(px, py) -(px.^2+py.^2); %-k^2
            obj.koperator = (bsxfun(kOper2, sample.grid.px_range, sample.grid.py_range));
        
            obj.grid = sample.grid;
            obj.roi  = sample.roi;
            obj.lambda  = options.lambda;
            obj.x_range = sample.grid.x_range(obj.roi{2});
            obj.x_range = obj.x_range - obj.x_range(1);
            obj.y_range = sample.grid.y_range(obj.roi{1});
            obj.y_range = obj.y_range - obj.y_range(1);
        end
        
        function E = exec(obj, sources)
            tic;
            %% preallocate the loop variables
            E_prev = zeros(obj.grid.N); % field at timestep = n-1
            E      = zeros(obj.grid.N);  % field at timestep = n
            
            %% initialize_sources_2d (scaling source size to grid size)
            %todo: respect sparsity
            source = zeros(obj.grid.N);
            source(obj.roi{1}, obj.roi{2}) = sources;
            
            %% Check whether gpu computation option is enabled
            if obj.gpuEnabled                
                E_prev = gpuArray(E_prev);
                E = gpuArray(E);
                source = gpuArray(source);
                obj.coeff = gpuArray(obj.coeff);
                obj.koperator = gpuArray(obj.koperator);
            end
                     
            
            %% Energy thresholds (convergence and divergence criterion)
            en_all    = zeros(1, obj.max_iterations);
            en_all(1) = wavesim.energy(source);
            threshold = obj.energy_threshold * en_all(1); %
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% PSTD time marching loop
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.it = 1;
            omega  = 2*pi/obj.lambda; %c==1
            while abs(en_all(obj.it)) >= threshold && obj.it <= obj.max_iterations
                obj.it = obj.it+1;
                
                %calculate source amplitude
                osc = obj.it * obj.dt * omega; 
                source_amplitude = exp(-1.0i * osc);
                %if (osc < 10*pi)
                %    source_amplitude = source_amplitude * (1-0.5*cos(osc/10));
                %end
                
                %update fields
                %based on equation:
                % nabla^2 E - e_r d^2 E/dt^2 = -source
                %
                %approximate time derivative by centered finite difference
                %
                % nabla^2 E - e_r (E_next - 2*E + E_prev) / dt^2 = -source
                %
                % isolate E_next:
                % E_next = (nabla^2 E + source) * dt^2/e_r + 2*E - E_prev
                E_next = obj.c2.*E + obj.c1.*E_prev + obj.c3 .* (ifft2(obj.koperator.*fft2(E)) + source_amplitude*source);
   
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%capture sampled wave fields_2d;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
                phase_shift = exp(-1.0i*obj.dt*omega); %expected phase shift for single step
                en_all(obj.it) = wavesim.energy(E_next-E*phase_shift);
                if (mod(obj.it, obj.callback_interval)==0) %now and then, call the callback function to give user feedback
                    obj.callback(obj, E, en_all(1:obj.it), threshold);
                end
                E_prev = E;
                E      = E_next;
            end %end timestep
            
            E = gather(E(obj.roi{1}, obj.roi{2}))/source_amplitude; % converts gpu array back to normal array and normalize phase shift from source
            obj.time = toc;
            disp(['Reached steady state in ' num2str(obj.it) ' iterations']);
            disp(['Time consumption: ' num2str(obj.time) ' s']);
        end
        
    end
    
    methods(Static)
        function options = readout_input(options)
           % function reading out all given constructor options and filling
           % in default values for missing properties
           % wavelength in real(V) = 0
           if ~isfield(options,'lambda')
               options.lambda = 1; % in um
           end
           % size of pixels
           if ~isfield(options,'pixel_size')
               options.pixel_size = 1/4*options.lambda; % in um
           end
        end
        
%         %default callback function. Shows real value of field
% 		function default_callback(obj, E, energy, threshold)
%             figure(1);
%             imagesc(obj.grid.x_range, obj.grid.y_range, abs(E));
%             title(['iteration = ' num2str(obj.it)])
%             xlabel('x (m)','FontSize',14); ylabel('y (m)','FontSize',14);
%             h = colorbar; set(get(h,'Title'),'String','Re(E) (a.u.)','FontSize',14);
%             set(gca,'FontSize',14);
%             drawnow;
%         end
        
        %default callback function. Shows real value of field, and total energy evolution
        function default_callback(obj, E, energy, threshold)
            figure(1);
            subplot(2,1,1); plot(1:length(energy),log10(energy),'b',[1,length(energy)],log10(threshold)*ones(1,2),'--r');
            title(length(energy));  xlabel('# iterations'); ylabel('log(energy added)');
            
            subplot(2,1,2); plot(real(E(end/2,:))); title('midline cross-section')
            xlabel('y (\lambda / 4)'); ylabel('real(E_x)');
            
            disp(['Added energy ', num2str(energy(end))]);
            drawnow;
        end
    end
    
end
    
    


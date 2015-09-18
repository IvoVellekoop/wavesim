classdef PSTD
    %Simulation of the 2-D wave equation using PSTD
    % Saroch Leedumrongwatthanakun 2015
    
    properties
        coeff;   % coefficient
        info; % diagnostics information on cover of field on medium
		gpuEnabled = false; % logical to check if simulation are ran on the GPU (default: false)
        callback = @PSTD.default_callback; %callback function that is called for showing the progress of the simulation. Default shows image of the real value of the field.
        callback_interval = 100; %the callback is called every 'callback_interval' steps. Default = 5
        dt; % time step
        time; % time
        number_of_time_steps; % a number of time step
    %%internal
        k; % wave number
        koperator; % k domain
        dtmin; %samll time step at sutible pixel_size
    end
    
    properties (Constant)
        eps_0 = 8.854187817e-12; % permittivity of free space               
        mu_0  = 4*pi*1e-7; % permeability of free space   
        c = 1/sqrt(PSTD.mu_0*PSTD.eps_0); % speed of light in free space
        c2 = 1/(PSTD.mu_0*PSTD.eps_0); % squre of speed of light in free space
    end
 
    methods
        function obj = PSTD(sample, options)
            %% Constructs a wave simulation object
			%	refractive_index = refractive index map (need not have an average value of 1)
			%	options.pixel_size = size of a pixel in the refractive_index map, e.g. in um
			%   options.wavelength = free space wavelength (same unit as pixel_size, e. g. um)
            %   options.time_size = size of time step in each loop in second
            
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use
            
            options = PSTD.readout_input(options);
                        
            %% Stability condition
            obj.dtmin = (2/sqrt(2)*sample.grid.dx/pi/PSTD.c); %PSTD
            %dt = 1/(c*sqrt((1/dx^2)+(1/dy^2))); %FDTD
            if (options.dt < obj.dtmin)
                obj.dt = options.dt;
            else
                obj.dt = 0.9*obj.dtmin; %add 0.9 factor to decrease dt
            end  
            
            % Time array
            obj.number_of_time_steps = round(options.time_duration/obj.dt)+2; %+2 for initial time and time at final step
            obj.time = (0:obj.number_of_time_steps+1)*obj.dt;
   
            %% Initialize updating coefficients;
            %obj.coeff = fftshift(obj.c2*options.dt^2./sample.refractive_index);
            obj.coeff = obj.c2*options.dt^2./(sample.refractive_index.^2);
            obj.coeff = obj.coeff(1:size(sample.damping_y,2),1:size(sample.damping_x,2)) .* (sample.damping_y' * sample.damping_x); 
            
            %%  Calculate PSTD operator
            %kOper1 = @(px, py) 1.0i*(px+py); %ik
            %obj.k = (bsxfun(kOper1, sample.grid.px_range, sample.grid.py_range));

            
            kOper2 = @(px, py) -(px.^2+py.^2); %-k^2
            obj.koperator = (bsxfun(kOper2, sample.grid.px_range, sample.grid.py_range));
            %obj.koperator = - (sample.grid.py_range * sample.grid.px_range).^2;
            
          
    
        end
        
        function [En, SumE, RMS, avg_time, n_time, success] = exec(sample,obj, source)
            %% preallocate the loop variables
            Enm = zeros(sample.grid.N(1),sample.grid.N(2)); % field at timestep = n-1
            En = zeros(sample.grid.N(1),sample.grid.N(2));  % field at timestep = n
            Enp = zeros(sample.grid.N(1),sample.grid.N(2)); % field at timestep = n+1
            SumE = zeros(sample.grid.N(1),sample.grid.N(2)); %sum of filed En
            RMS = zeros(sample.grid.N(1),sample.grid.N(2));
            
            
            %% initialize_sources_2d (scaling source size to grid size)
            source.value = padarray(source.value,sample.grid.N - size(source.value), 'post');
            %source = (1/PSTD.eps_0)*ifft2(bsxfun(@times, obj.k, fft2(source)));
            %source = circshift(source,round(-(2*boundaries.width+sample.grid.padding)/2));
             
            % copy waveform of the waveform type to waveform of the source 
            source.timeseries = source.magnitude_factor * source.timeseries; 
            %source.timeseries = source.magnitude_factor * source.WaveForm(source,obj.time); 
            %source.waveform = source.magnitude_factor * waveforms.source.waveform_type.waveform;
            
            %% Check whether gpu computation option is enabled
            if obj.gpuEnabled                
                Enm = gpuArray(Enm);
                En = gpuArray(En);
                Enp = gpuArray(Enp);
                SumE = gpuArray(SumE);
                RMS = gpuArray(SumE);
                source = gpuArray(source);
                obj.coeff = gpuArray(obj.coeff);
                obj.koperator = gpuArray(obj.koperator);
                %Enm = source.value;
                %En = source.value;
            else
                %Enm = source.value;
                %En = source.value;
            end
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% PSTD time marching loop
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            success = true;
            start_time = cputime; 
            current_time = 0;
            avg_time = 0;
            cover_time = max(sample.N)*sample.grid.dx/PSTD.c; % time of wave propagation thorug whole medium
             %source.timeseries(round(cover_time/obj.dt/4):end) = 0; 
            n_time=0;
            for time_step = 2:obj.number_of_time_steps-1
                
                current_time  = current_time + obj.dt;
                
                %%update_source;
                %source_term =0;
                %source_term = source.timeseries(time_step)*source.value;  
                source_term = (source.timeseries(time_step+1) - source.timeseries(time_step-1))/2/obj.dt*source.value; 
                %En = En + source.timeseries(time_step)*source.value; 
           
                %%update_fields_2d;
                Enp = 2*En - Enm ...
                        + bsxfun(@times, obj.coeff , ...
                            (-source_term + (ifft2( ...
                                              bsxfun(@times, obj.koperator, fft2( En ))))));
                                          
                %Enp = 2*En - Enm + obj.coeff .*(-source_term + (ifft2( obj.koperator.*fft2(En) )));
                
   
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%capture sampled wave fields_2d;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
               	if (mod(time_step, obj.callback_interval)==0) %now and then, call the callback function to give user feedback
                    obj.callback(En,sample,current_time);
                    disp(['time_step ', num2str(time_step)]);
                end
                if (current_time > cover_time)
                    RMS = RMS + (En.^2*obj.dt);
                    SumE = SumE + real(En*exp(2.0i*pi*source.sinusoidal.frequency*current_time));
                    n_time= n_time+1;
                    avg_time  = avg_time + obj.dt;
                end
                
                %% updating field
                %En(1:sample.N(1),1) = 1;
                %Enp(1:sample.N(1),1) = 1;
                Enm = En;
                En = Enp;
    
            end %end timestep
            SumE = SumE / n_time;
            RMS = sqrt ( RMS / avg_time);
            end_time = cputime;
            total_time_in_minutes = (end_time - start_time)/60;
            disp(['Total simulation time is ' num2str(total_time_in_minutes) ' minutes.']);
            
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
        
        %default callback function. Shows real value of field
		function default_callback(E,sample,current_time)
            figure(1);
            imagesc(sample.grid.x_range,sample.grid.y_range,real(E(1:sample.N(1),2:sample.N(2))));
            title(['time =' num2str(current_time) ' s'])
            xlabel('x (m)','FontSize',14); ylabel('y (m)','FontSize',14);
            h = colorbar; set(get(h,'Title'),'String','Re(E) (a.u.)','FontSize',14);
            set(gca,'FontSize',14);
            drawnow;
        end
    end
    
end
    
    


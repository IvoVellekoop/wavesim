classdef wavesim
    %Simulation of the 2-D wave equation using a Born series approach
    % Ivo M. Vellekoop 2014
    
    properties
        V;   % potential array used in simulation
        info; % diagnostics information on convergence
		gpuEnabled = false; % logical to check if simulation are ran on the GPU (default: false)
        callback = @wavesim.default_callback; %callback function that is called for showing the progress of the simulation. Default shows image of the absolute value of the field.
		callback_interval = 1000; %the callback is called every 'callback_interval' steps. Default = 5
		energy_threshold = 1E-9; %the simulation is terminated when the added energy between two iterations is lower than 'energy_threshold'. Default 1E-9
		max_iterations = 1E4; %1E4; %or when 'max_iterations' is reached. Default 10000
        it; %iteration
        time; % time comsumption
     %%internal
        k; % wave number
        epsilon; 
        epsilonmin; 
        g0_k; % bare Green's function used in simulation
        %performance

    end
        
    methods
        function obj=wavesim(sample, options)
			%% Constructs a wave simulation object
			%	refractive_index = refractive index map (need not have an average value of 1)
			%	options.pixel_size = size of a pixel in the refractive_index map, e.g. in um
			%   options.wavelength = free space wavelength (same unit as pixel_size, e. g. um)
            %   options.epsilon = convergence parameter

                        
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use
            %% Read out constructor options and fill out missing properties with defualt values
            if nargin == 1
                sample = struct;
                options = struct;
            end
            options = wavesim.readout_input(options);
            
            
            %% Determine constants based on refractive_index map
            obj.k = sample.n_center * (2*pi/options.lambda);   
            
            %% Stability condition
            fluctuation = abs((sample.refractive_index(1:sample.N(1),1:sample.N(2))*2*pi/options.lambda).^2-obj.k^2).^2;
            % old version:(2*pi/options.lambda)^2 * (n_max^2 - n_min^2)/2
            obj.epsilonmin = max(sqrt(fluctuation(:)));
            if (options.epsilon > obj.epsilonmin)
                obj.epsilon = options.epsilon;
            else
                obj.epsilon = max(sample.BCepsmin,obj.epsilonmin); %compare epsilon from BC and MoBorn
            end  
            
            %% Calculate Green function for k_red (reduced k vector: k_red^2 = k_0^2 + 1.0i*epsilon)
            f_g0_k = @(px, py) 1./(px.^2+py.^2-(obj.k^2 + 1.0i*obj.epsilon));
            obj.g0_k = bsxfun(f_g0_k, sample.grid.px_range, sample.grid.py_range);
            
            %% Potential map (V==k^2-k_0^2-1i*epsilon). (First pad refractive index map)            
            obj.V = (sample.refractive_index*2*pi/options.lambda).^2-(obj.k^2 + 1.0i*obj.epsilon);
            
            %% Apply absorbing boundaries on potential map
            obj.V = obj.V(1:size(sample.damping_y,2),1:size(sample.damping_x,2)) .* (sample.damping_y' * sample.damping_x); 
               
        end
        
        function epsilon = epsiloncondition(obj)
           disp(['the minimum value of epsilon is ' num2str(obj.epsilonmin)]);
        end
               
        function [E_x, iter, time, success] = exec(sample,obj, source)
            %%% Execute simulation 
            time=0;
            iter=1;
            tic;
            %% Scaling source size to grid size
            source = padarray(source,sample.grid.N - size(source), 'post');

            %% Check whether gpu computation option is enabled
            if obj.gpuEnabled                
                obj.g0_k = gpuArray(obj.g0_k);
                obj.V = gpuArray(obj.V);
                source = gpuArray(source);
                en_all = zeros(1, obj.max_iterations,'gpuArray');
                E_x = ifft2(obj.g0_k .* fft2(source));
            else
                E_x = ifft2(obj.g0_k .* fft2(source));
                en_all = zeros(1, obj.max_iterations);
            end
  
            %% Energy thresholds (convergence and divergence criterion)
            en_all(1) = wavesim.energy(E_x);
            threshold = obj.energy_threshold * en_all(1); % 
            div_thresh = 10 * en_all(1);
            
            %% simulation iterations
            success = true;
            obj.it = 1;            
            while abs(en_all(obj.it)) >= threshold && obj.it <= obj.max_iterations
                obj.it = obj.it+1;
                Eold = E_x;
                E_x = single_step(obj, E_x, source);
                en_all(obj.it) = wavesim.energy(E_x - Eold);
				if (mod(obj.it, obj.callback_interval)==0) %now and then, call the callback function to give user feedback
                    obj.callback(E_x(end/2,:), en_all(1:obj.it), threshold);
                end
                
                if abs(en_all(obj.it)) > div_thresh || isnan(en_all(obj.it)) % abort when added energy is more than 10 times the initial added energy
                    success = false;
                    break;
                end                               
            end
            %E_x = (1.0i*obj.V/obj.epsilon).*gather(E_x); % converts gpu array back to normal array
            E_x = gather(E_x); % converts gpu array back to normal array
            %% Simulation finished
            obj.time=toc;
            time=obj.time;
            iter=obj.it;
            if success && abs(en_all(obj.it)) < threshold
                disp(['Reached steady state in ' num2str(obj.it) ' iterations']);
                disp(['Time comsumption: ' num2str(obj.time) ' s']);
            elseif ~success
                disp('simulation diverged');
            end
        end

        function [E_x, energy_diff] = single_step(obj, E_x, source) 
            % performs a single iteration of the algorithm
            E_x = E_x - (1.0i*obj.V/obj.epsilon) .* (E_x-ifft2(obj.g0_k .* fft2(obj.V.*E_x+source))); %wavesim version 
            %E_x = E_x - (1.0i*obj.V/obj.epsilon) .* (E_x-ifft2(obj.g0_k .*
            %fft2(obj.V.*E_x)))+ifft2(obj.g0_k.*fft2(source)); %cutout preconditioner on source
            %E_x = E_x - (1.0i*obj.V/obj.epsilon).*E_x+ifft2(obj.g0_k.*(fft2(obj.V.*(1.0i*obj.V/obj.epsilon).* E_x))) + ifft2(obj.g0_k .* fft2(source)); %version 6
        end
        
        function [E_propagation, en_all , time, success] = psudopropagation(obj, source, numberiter,frame)
            %%% Execute simulation 
            time=0;
            tic;
            %% Scaling source size to grid size
            source = padarray(source,obj.grid.N - size(source), 'post');

            %% Check whether gpu computation option is enabled
            if obj.gpuEnabled                
                obj.g0_k = gpuArray(obj.g0_k);
                obj.V = gpuArray(obj.V);
                source = gpuArray(source);
                en_all = zeros(1, numberiter+1,'gpuArray');
                E_x = ifft2(obj.g0_k .* fft2(source));
                %% Propagation of field parameter
                E_propagation=gpuArray(zeros(size(E_x,1),size(E_x,2),round(numberiter/frame)));
                E_propagation(:,:,1)=E_x; %inital field
            else
                %{
                %filtering high frequency of source term
                source_k=fft2(source);
                source_k(1,3*end/4:end)=0;
                source_k(1,1:end/4)=0;
                source = ifft2(source_k);
                %}
                E_x = ifft2(obj.g0_k .* fft2(source));
                en_all = zeros(1, numberiter+1);
                %% Propagation of field parameter
                E_propagation=zeros(size(E_x,1),size(E_x,2),round(numberiter/frame));
                E_propagation(:,:,1)=E_x; %inital field
            end
        
            %% Energy thresholds (convergence and divergence criterion)
            en_all(1) = wavesim.energy(E_x);
            threshold = obj.energy_threshold * en_all(1); % 
            div_thresh = 10 * en_all(1);
            
            %% simulation iterations
            success = true;
            obj.it = 1;            
            while en_all(obj.it) >= threshold && obj.it <= numberiter
                obj.it = obj.it+1;
                [E_x, diff_energy] = single_step(obj, E_x, source);
                
                %E_propagation(obj.it,:)=E_x;
                en_all(obj.it) = diff_energy;
                
                if (mod(obj.it, frame)==0) %now and then, call the callback function to give user feedback
                    E_propagation(:,:,(obj.it/frame))=E_x;
                end;
				if (mod(obj.it, obj.callback_interval)==0) %now and then, call the callback function to give user feedback
                    disp(['Iteration loop:' num2str(obj.it)]);
                    obj.callback(E_x(end/2,:), en_all(1:obj.it), threshold);
				end;
                
                if abs(en_all(obj.it)) > div_thresh || isnan(en_all(obj.it)) % abort when added energy is more than 10 times the initial added energy
                    success = false;
                    break;
                end
            end;
            E_propagation=gather(E_propagation);  % converts gpu array back to normal array
            %% Simulation finished
            obj.time=toc;
            time=obj.time;
            if success && abs(en_all(obj.it)) < threshold
                disp(['Reached steady state in ' num2str(obj.it) ' iterations']);
                disp(['Time comsumption: ' num2str(obj.time) ' s']);
            elseif ~success
                disp('simulation diverged');
            end
        end
                       
        function analyze(obj)
            %% Displays various information
            g = obj.grid;
            inf = obj.info;
            disp(['Size of simulation: ', num2str(g.N*g.dx), ' wavelengths']); %num2str(g.Nx*g.dx)
            disp(['Number of samples per wavelength: ', num2str(1/g.dx)]);
            disp(['Relative bandwidth reserved for refractive index map: ', num2str(obj.bandwidth)]); % unknown bandwidth
            disp(['Smallest feature in refractive index map: ', num2str(g.dx/obj.bandwidth), ' wavelengths']); % unknown bandwidth
            disp(['Smallest feature in field: ', num2str(g.dx/(1-obj.bandwidth)), ' wavelengths']); % unknown bandwidth
            disp(['Truncation of Green function causes loss factor of ', num2str(1-min(inf.truncation_loss_P, inf.truncation_loss_g0_k_max)), ' per step']);

            %% Given an estimate for the stability of the simulation
            %% estimate exponent of divergence by calculating local average of V
            % the local average is calculated by convolving with |g_x|
            g_x_abs_f = fft2(abs(ifft2(obj.g0_k)));
            Vconv = ifft2(fft2(obj.V).*g_x_abs_f/g_x_abs_f(1,1));
            Vconv_max = max(abs(Vconv(:)));
            %divergence_homogeneous = abs(p_factor * g0_k_max); 
            divergence_worst_case = abs(Vconv_max * inf.trunc_g0_k_max); 
            %disp(['Divergence exponent in homogeneous medium ', num2str(divergence_homogeneous)]);
            disp(['Divergence exponent, worst case ', num2str(divergence_worst_case)]);
            disp(['Total divergence single pass ', num2str(divergence_worst_case^(min(g.Nx, g.Ny)*g.dx/obj.epsilon),'%E')]);            
        end
        %callback function for performance
		function [it time] = performance(obj)
                it=obj.it;
                time=obj.time;
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
        
        function en = energy(E_x)
            en= sum(abs(E_x(:)).^2);
        end
		
		%default callback function. Shows real value of field, and total energy evolution
		function default_callback(E_cross, energy, threshold)
            figure(1);
			subplot(2,1,1); plot(1:length(energy),log10(energy),'b',[1,length(energy)],log10(threshold)*ones(1,2),'--r'); 
            title(length(energy));  xlabel('# iterations'); ylabel('log(energy added)');
            
            subplot(2,1,2); plot(real(E_cross)); title('midline cross-section')
            xlabel('y (\lambda / 4)'); ylabel('real(E_x)');
            			
            disp(['Added energy ', num2str(energy(end))]); 
			drawnow;
        end
    end
end
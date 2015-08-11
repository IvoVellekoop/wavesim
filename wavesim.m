classdef wavesim
    %Simulation of the 2-D wave equation using a Born series approach
    % Ivo M. Vellekoop 2014
    
    properties
        V;   % potential array used in simulation
        info; % diagnostics information on convergence
        grid; % simgrid object
		gpuEnabled = false; % logical to check if simulation are ran on the GPU (default: false)
        callback = @wavesim.default_callback; %callback function that is called for showing the progress of the simulation. Default shows image of the absolute value of the field.
		callback_interval = 100; %the callback is called every 'callback_interval' steps. Default = 5
		energy_threshold = 1E-4; %the simulation is terminated when the added energy between two iterations is lower than 'energy_threshold'. Default 1E-9
		max_iterations = 1E4; %1E4; %or when 'max_iterations' is reached. Default 10000
        it; %iteration
        time; % time comsumption
     %%internal
        k; % wave number
        epsilon; % 1/'window size'
        g0_k; % bare Green's function used in simulation
        %performance

    end
        
    methods
        function obj=wavesim(refractive_index, options, boundaries)
			%% Constructs a wave simulation object
			%	refractive_index = refractive index map (need not have an average value of 1)
			%	options.pixel_size = size of a pixel in the refractive_index map, e.g. in um
			%   options.wavelength = free space wavelength (same unit as pixel_size, e. g. um)
            %   options.epsilon = convergence parameter
            %	boundaries.width = extra space to add around simulation to simulate absorbing boundaries, in pixels
			%   boundaries.xcurve = shape of the damping curve horizontal boundary
            %   boundaries.ycurve = shape of the damping curve vertical boundary
                        
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use
            %% Read out constructor options and fill out missing properties with defualt values
            if nargin == 1
                options = struct;
                boundaries = struct;
            elseif nargin == 2
                boundaries = struct;
            end
            [options,boundaries] = wavesim.readout_input(options,boundaries,size(refractive_index));
            
            %% setup grid, taking into account required boundary. Pad to next power of 2 when needed
            obj.grid = simgrid(size(refractive_index)+2*boundaries.width, options.pixel_size);
            
            %% Padding refractive index to match simulation grid           
            refractive_index = padarray(refractive_index, round((2*boundaries.width+obj.grid.padding)/2), 'replicate', 'both');
            refractive_index = circshift(refractive_index,round(-(2*boundaries.width+obj.grid.padding)/2));

            %% Determine constants based on refractive_index map
			n_min = min(abs(refractive_index(:))); 
            n_max = max(abs(refractive_index(:))); 
            n_center = sqrt((n_max^2 + n_min^2) / 2); %central refractive index (refractive index that k_r is based on)
            obj.k = n_center * (2*pi/options.lambda);
            
			%% Determine optimum value for epsilon, otherwise use given epsilon (epsilon = 1/step size)
            if isfield(options,'epsilon') % check if epsilon is given as input
                obj.epsilon = options.epsilon;
            else
                epsmin = (2*pi/options.lambda) / (max(boundaries.width)*options.pixel_size); %epsilon cannot be smaller, or green function would wrap around boundary (pre-factor needed!)
                obj.epsilon = 1.1*max(epsmin, (2*pi/options.lambda)^2 * (n_max^2 - n_min^2)/2); %the factor 1.1 is a safety margin to account for rounding errors.
            end  
    
			%% Calculate Green function for k_red (reduced k vector: k_red^2 = k_0^2 + 1.0i*epsilon)
            f_g0_k = @(px, py) 1./(px.^2+py.^2-(obj.k^2 + 1.0i*obj.epsilon));
            obj.g0_k = bsxfun(f_g0_k, obj.grid.px_range, obj.grid.py_range);
            
            %% Potential map (V==k^2-k_0^2-1i*epsilon). (First pad refractive index map)            
            obj.V = (refractive_index*2*pi/options.lambda).^2-(obj.k^2 + 1.0i*obj.epsilon);                        
            
            %% Absorbing boundaries
            damping_x = [ ones(1, obj.grid.N(2)-obj.grid.padding(2)-2*boundaries.width(2)), boundaries.xcurve, zeros(1, obj.grid.padding(2)), fliplr(boundaries.xcurve)];
            damping_y = [ ones(1, obj.grid.N(1)-obj.grid.padding(1)-2*boundaries.width(1)), boundaries.ycurve, zeros(1, obj.grid.padding(1)), fliplr(boundaries.ycurve)];
                        
            %% Apply absorbing boundaries on potential map
            obj.V = obj.V(1:size(damping_y,2),1:size(damping_x,2)) .* (damping_y' * damping_x);
            
        end;       
               
        function [E_x, iter, time, success] = exec(obj, source)
            %%% Execute simulation 
            time=0;
            iter=1;
            tic;
            %% Scaling source size to grid size
            source = padarray(source,obj.grid.N - size(source), 'post');

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
                [E_x, diff_energy] = single_step(obj, E_x, source);
                en_all(obj.it) = diff_energy;
                
				if (mod(obj.it, obj.callback_interval)==0) %now and then, call the callback function to give user feedback
                    obj.callback(E_x(:,end/2), en_all(1:obj.it), threshold);
				end;
                
                if abs(en_all(obj.it)) > div_thresh || isnan(en_all(obj.it)) % abort when added energy is more than 10 times the initial added energy
                    success = false;
                    break;
                end                               
            end;
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
            % returns difference energy (todo: optimize difference energy
            % calculation)
            Eold = E_x;
            E_x = E_x - (1.0i*obj.V/obj.epsilon) .* (E_x-ifft2(obj.g0_k .* fft2(obj.V.*E_x+source))); %wavesim version 
            %E_x = E_x - (1.0i*obj.V/obj.epsilon) .* (E_x-ifft2(obj.g0_k .*
            %fft2(obj.V.*E_x)))+ifft2(obj.g0_k.*fft2(source)); %cutout preconditioner on source
            %E_x = E_x - (1.0i*obj.V/obj.epsilon).*E_x+ifft2(obj.g0_k.*(fft2(obj.V.*(1.0i*obj.V/obj.epsilon).* E_x))) + ifft2(obj.g0_k .* fft2(source)); %version 6
            energy_diff = wavesim.energy(E_x-Eold);
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
                    obj.callback(E_x(:,end/2), en_all(1:obj.it), threshold);
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
        function [options,boundaries] = readout_input(options,boundaries,size)
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
           
           % fill in boundary settings
           if ~isfield(boundaries,'width') % if no bonudary parameters are given
               boundaries.width = floor( (2.^nextpow2(size) - size)/4); % width of the absorbing boundary
               boundaries.xcurve = 1-linspace(0, 1, boundaries.width(2)).^2;   % damping curve of horizontal absorbing boundary
               boundaries.xcurve = 1-linspace(0, 1, boundaries.width(2)).^2;   % curve of vertical absorbing layer
           end
                    
           if ~isfield(boundaries,'xcurve') || boundaries.width(2) ~= length(boundaries.xcurve) % if width is given but not xcurve
               boundaries.xcurve = 1-linspace(0, 1, boundaries.width(2)).^2;
           end
           
           if ~isfield(boundaries,'ycurve') || boundaries.width(1) ~= length(boundaries.ycurve) % if width is given but not ycurve
               boundaries.ycurve = 1-linspace(0, 1, boundaries.width(1)).^2;
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
        end;
    end
end
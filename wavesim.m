classdef wavesim
    %Simulation of the 2-D wave equation using a Born series approach
    % Ivo M. Vellekoop 2014
    
    properties
        V;   % potential array used in simulation
        info; % diagnostics information on convergence
        grid; % simgrid object
        bandwidth;
        background; % background mask
		gpuEnabled; % logical to check if simulation are ran on the GPU (default: false)
        callback; %callback function that is called for showing the progress of the simulation. Default shows image of the absolute value of the field.
		callback_interval; %the callback is called every 'callback_interval' steps. Default = 5
		energy_threshold; %the simulation is terminated when the added energy between two iterations is lower than 'energy_threshold'. Default 1E-9
		max_iterations; %or when 'max_iterations' is reached. Default 10000
     %%internal
        k; % wave number
        k_red2; % k0^2+i*epsilon
        epsilon; % 1/'window size'
        g0_k; % bare Green's function used in simulation
    end
        
    methods
        function obj=wavesim(refractive_index, pixel_size, lambda, boundary, bandwidth, epsilon)
			%% Constructs a wave simulation object
			%	refractive_index = refractive index map (need not have an average value of 1)
			%	pixel_size = size of a pixel in the refractive_index map, e.g. in um
			%   wavelength = free space wavelength (same unit as pixel_size, e. g. um)
			%	boundary = extra space to add around simulation to simulate absorbing boundaries, in pixels 
			%
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use
			
			%% set default options
			obj.callback = @wavesim.default_callback;
			obj.callback_interval = 1000;
            obj.energy_threshold = 1E-9; % fraction of initially added energy
			obj.max_iterations = 200000;
            obj.gpuEnabled = false;
            
            %% Determine constants based on refractive_index map
			n_min = min(refractive_index(:));
            n_max = max(refractive_index(:));
            n_center = sqrt((n_max^2 + n_min^2) / 2); %central refractive index (refractive index that k_r is based on)
            obj.k = n_center * (2*pi/lambda);
            
			%% determine optimum value for epsilon (epsilon = 1/step size)
            epsmin = (2*pi/lambda) / (boundary*pixel_size); %epsilon cannot be smaller, or green function would wrap around boundary (pre-factor needed!)
            obj.epsilon = max(epsmin, (2*pi/lambda)^2 * (n_max^2 - n_min^2)/2) * 1.3; %the factor 1.1 is a safety margin to account for rounding errors.
            
            if exist('epsilon','var') % check if epsilon is given
                obj.epsilon = epsilon;
            end
            
            obj.k_red2 = obj.k^2 + 1.0i*obj.epsilon; %k reduced squared
                        
			%% setup grid, taking into account required boundary. Pad to next power of 2 when needed
			obj.grid = simgrid(size(refractive_index)+2*boundary, pixel_size);
            
			%% Calculate Green function for k_red (reduced k vector)
            f_g0_k = @(px, py) 1./(px.^2+py.^2-obj.k_red2);
            obj.g0_k = bsxfun(f_g0_k, obj.grid.px_range, obj.grid.py_range);
            
            %% Potential map (V==k_r^2-k^2). (First pad refractive index map)
            refractive_index = padarray(refractive_index, obj.grid.N-size(refractive_index), n_center, 'post');
            obj.V = (refractive_index*2*pi/lambda).^2-obj.k_red2;
             
            %% Low pass filter potential function V
            obj.bandwidth = bandwidth;
            width = round(min(obj.grid.N)*bandwidth/2)*2;
            win2d = tukeywin(width, 0.125) * tukeywin(width, 0.125).';
            win2d = fftshift(padarray(win2d, (obj.grid.N-size(win2d))/2, 0));
            obj.V = ifft2(win2d.*fft2(obj.V));
            
			%% defining damping curve
            d_curve = 1-linspace(0, 1, boundary).^2;
            damping_x = [ ones(1, obj.grid.N(2)-obj.grid.padding(2)-2*boundary), d_curve, zeros(1, obj.grid.padding(2)), d_curve(end:-1:1)];
            damping_y = [ ones(1, obj.grid.N(1)-obj.grid.padding(1)-2*boundary), d_curve, zeros(1, obj.grid.padding(1)), d_curve(end:-1:1)];
            
            obj.V = obj.V .* (damping_y' * damping_x);
            obj.background = damping_y' * damping_x < 1;
        end;
            
        function [E_x,success] = exec(obj, source)
            %% Execute simulation 
            % Scaling source size to grid size
            [a,b] = size(source);
            source = padarray(source,obj.grid.N - [a,b], 'post');
            source(obj.background) = 0;

            %% Check whether gpu computation option is enabled
            if obj.gpuEnabled
                E_x = zeros(size(source),'gpuArray');
                E_k = zeros(size(source),'gpuArray');
                obj.g0_k = gpuArray(obj.g0_k);
                obj.V = gpuArray(obj.V);
                source = gpuArray(full(source)); % gpu does not support sparse matrices
                en_all = zeros(1, obj.max_iterations,'gpuArray');
            else
                E_x = 0;
                en_all = zeros(1, obj.max_iterations);
            end
            
            success = true;
            %% simulation iterations
             for it=1:obj.max_iterations
                E_old = E_x; % previous iteration
                E_ak = obj.g0_k .* fft2(obj.V.*E_x + source);
                E_x = ifft2(E_ak) + conj(obj.V) .* ifft2(conj(obj.g0_k) .* (fft2(E_x) - E_ak));
                
                en_all(it) = wavesim.energy(E_x - E_old);
                
                if it == 1 % after first iteration the energy threshold is determined
                    threshold = obj.energy_threshold * en_all(1);
                end
                
				if (mod(it, obj.callback_interval)==0) %now and then, call the callback function to give user feedback					
                    obj.callback(E_x(:,round(b/2)), en_all(1:it), threshold);
				end;
                
                if it>100 && (abs(en_all(it) > 10 * en_all(1) || isnan(abs(en_all(it))))) % abort when added energy is more than 10 times the initial added energy
                    disp('simulation diverged');
                    success = false;
                    break;
                end
                    
                if it>100 && abs(en_all(it)) < obj.energy_threshold * en_all(1) % simulation finished when energy threshold is reached
                    disp(['Reached steady state in ' num2str(it) ' iterations']);
                    break;
                end;
            end;
            E_x = gather(E_x);
        end;
        
        function analyze(obj)
            %% Displays various information
            g = obj.grid;
            inf = obj.info;
            disp(['Size of simulation: ', num2str(g.Nx*g.dx), ' wavelengths']);
            disp(['Number of samples per wavelength: ', num2str(1/g.dx)]);
            disp(['Relative bandwidth reserved for refractive index map: ', num2str(obj.bandwidth)]);
            disp(['Smallest feature in refractive index map: ', num2str(g.dx/obj.bandwidth), ' wavelengths']);
            disp(['Smallest feature in field: ', num2str(g.dx/(1-obj.bandwidth)), ' wavelengths']);
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
        end;
    end
    methods(Static)
        function en = energy(E_x)
            en= sum(abs(E_x(:)).^2);
        end;
		
		%default callback function. Shows real value of field, and total energy evolution
		function default_callback(E_cross, energy, threshold)		    
			subplot(2,1,1); plot(1:length(energy),log10(energy),'b',[1,length(energy)],log10(threshold)*ones(1,2),'--r'); 
            title(length(energy));  xlabel('# iterations'); ylabel('log(energy added)');
            
            subplot(2,1,2); plot(real(E_cross)); title('midline cross-section')
            xlabel('y (\lambda / 4)'); ylabel('real(E_x)');
			
            disp(['Added energy ', num2str(energy(end))]); 
			drawnow;
        end;
    end
end
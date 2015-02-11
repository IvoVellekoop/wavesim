classdef wavesim
    %Simulation of the 2-D wave equation using a Born series approach
    % Ivo M. Vellekoop 2014
    
    properties
        V;   % potential array used in simulation
        info; % diagnostics information on convergence
        grid; % simgrid object
        bandwidth;
        background; % background mask
        g_curve;
		callback; %callback function that is called for showing the progress of the simulation. Default shows image of the absolute value of the field.
		callback_interval; %the callback is called every 'callback_interval' steps. Default = 5
		energy_threshold; %the simulation is terminated when the added energy between two iterations is lower than 'energy_threshold'. Default 1E-9
		max_iterations; %or when 'max_iterations' is reached. Default 10000
     %%internal
        k; % wave number
        k_red; % k-i*epsilon
        epsilon; % 1/'window size'
        g0_k; % bare Green's function used in simulation
    end
        
    methods
        function obj=wavesim(refractive_index, pixel_size, lambda, boundary, bandwidth)
			%% Constructs a wave simulation object
			%	refractive_index = refractive index map (need not have an average value of 1)
			%	pixel_size = size of a pixel in the refractive_index map, e.g. in um
			%   wavelength = free space wavelength (same unit as pixel_size, e. g. um)
			%	boundary = extra space to add around simulation to simulate absorbing boundaries, in pixels 
			%
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use
			
			%% set default options
			obj.callback = wavesim.default_callback;
			obj.callback_interval = 5;
            obj.energy_threshold = 1E-12;
			obj.max_iterations = 10000;

            %% Determine constants based on refractive_index map
			n_min = min(refractive_index(:));
            n_max = max(refractive_index(:));
            n_center = sqrt((n_max^2 + n_min^2) / 2); %central refractive index (refractive index that k_r is based on)
            obj.k = n_center * (2*pi/lambda);
            
			%% determine optimum value for epsilon (epsilon = 1/step size)
            epsmin = (2*pi/lambda) / (boundary*dx); %epsilon cannot be smaller, or green function would wrap around boundary (pre-factor needed!)
            obj.epsilon = max(epsmin, (2*pi/lambda)^2 * (n_max^2 - n_min^2)/2) * 1.1; %the factor 1.1 is a safety margin to account for rounding errors.
            obj.k_red = sqrt(obj.k^2 + 1.0i*obj.epsilon);
		
			%% setup grid, taking into account required boundary. Pad to next power of 2 when needed
			obj.grid = simgrid(n_size+2*boundary, pixel_size);
            
			%% Calculate Green function for k_red (reduced k vector)
            f_g0_k = @(px, py) 1./(px.^2+py.^2-obj.k_red^2);
            obj.g0_k = bsxfun(f_g0_k, obj.grid.px_range, obj.grid.py_range);
      
            %% Potential map (V==k_r^2-k^2). (First pad refractive index map)
            refractive_index = padarray(refractive_index, obj.grid.N-size(refractive_index), n_center, 'post');
            obj.V = ones(obj.grid.Ny, obj.grid.Nx) * (obj.k^2-obj.k_red^2);
            obj.V(1:size(refractive_index,1), 1:size(refractive_index,2)) = (refractive_index*2*pi/lambda).^2-obj.k_red^2;
            
             
            %% Low pass filter potential function V
            obj.bandwidth = bandwidth;
            width = round(min(obj.grid.N)*bandwidth);
            w_x = [zeros(ceil((obj.grid.N(2)-width)/2),1); tukeywin(width, 0.125); zeros(floor((obj.grid.N(2)-width)/2),1)];
            w_y = [zeros(ceil((obj.grid.N(1)-width)/2),1); tukeywin(width, 0.125); zeros(floor((obj.grid.N(1)-width)/2),1)].';
            win2d = fftshift(bsxfun(@times, w_x, w_y));
            obj.V = ifft2(win2d.*fft2(obj.V));
            
			%% defining damping curve
            obj.g_curve = 1-linspace(0, 1, obj.grid.x_padding).^2;
            damping_x = [ ones(1, obj.grid.N(2)-obj.grid.padding(2)), obj.g_curve, obj.g_curve(end:-1:1)];
            damping_y = [ ones(1, obj.grid.N(1)-obj.grid.padding(1)), obj.g_curve, obj.g_curve(end:-1:1)];
            obj.V = obj.V * damping_y' * damping_x;         
        end;
            
        function [E_x] = exec(obj, source)
            %% Execute simulation
            E_x = 0;
            source(obj.background) = 0;
                        
            en_all = zeros(1, obj.max_iterations);
            for it=1:obj.max_iterations
				%% perform the iteration
				%old method: E_x = ifft2(obj.g0_k .* fft2(obj.V.*E_x+source));
                E_a = ifft2(obj.g0_k .* fft2(obj.V.*E_x+source));
				E_x = ifft2(conj(obj.g0_k) .* fft2(conj(obj.V).*(E_x-E_a))) + E_a;
        
				en_all(it) = wavesim.energy(E_x);
				if (mod(it, obj.callback_interval)==0) %now and then, call the callback function to give user feedback
					obj.callback(E_x, en_all(1:it));
				end;
                if (abs(en_all(end)-en(all(end-1))) < obj.energy_threshold) %abort when threshold is reached
                    disp(['Reached steady state in ' num2str(it) ' iterations']);
                    break;
                end;
            end;            
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
            en= E_x(:)' * E_x(:);
        end;
		
		%default callback function. Shows real value of field, and total energy evolution
		function default_callback(E_x, energy)
		    subplot(2,1,2); plot(real(E_x(:,end/2)));
			subplot(2,1,1); plot(energy); title(it);
			disp(['Added energy ', num2str(energy(end)-energy(end-1))]); 
			drawnow;
        end;
    end
end
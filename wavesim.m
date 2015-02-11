classdef wavesim
    %Simulation of the 2-D wave equation using a Born series approach
    % Ivo M. Vellekoop 2014
    
    properties
        V;   % potential array used in simulation
        info; % diagnostics information on convergence
        grid; % simgrid object
        bandwidth;
        background; % background mask
        g_curve
     %%internal
        k; % wave number
        k_red; % k-i*epsilon
        epsilon; % 1/'window size'
        g0_k; % bare Green's function used in simulation
    end
        
    methods
        function obj=wavesim(refractive_index, bandwidth)
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use

            %% Determine constants based on refractive_index map
			dx = 1/4; %hard coded: 4 grid points per wavelength
            lambda = 1;
            n_min = min(refractive_index(:));
            n_max = max(refractive_index(:));
            n_center = sqrt((n_max^2 + n_min^2) / 2); %central refractive index (refractive index that k_r is based on)
            obj.k = n_center * (2*pi/lambda);
            
			%determine optimum value for epsilon (epsilon = 1/step size)
            epsmin = (2*pi/lambda)^2 * 0.05;
            obj.epsilon = max(epsmin, (2*pi/lambda)^2 * (n_max^2 - n_min^2)/2);
            obj.k_red = sqrt(obj.k^2 + 1.0i*obj.epsilon);
            
			%% make sure that size is a multiple of 2 (for simgrid), then setup coordinates
			n_size = ceil(size(refractive_index)/2)*2;
			refractive_index = padarray(refractive_index, size(refractive_index)-n_size, n_center, 'post');
            obj.grid = simgrid(n_size(2), n_size(1), dx);
            disp(['epsilon = ',num2str(obj.epsilon)]);
      
			%% Calculate Green function for k_red (reduced k vector)
            f_g0_k = @(px, py) 1./(px.^2+py.^2-obj.k_red^2);
            obj.g0_k = bsxfun(f_g0_k, obj.grid.px_range, obj.grid.py_range);

            %% truncate Green function in x-space to have a finite support
            % for the convolution.
            % This operation is slightly lossy!
            % First find the distance at which g_x drops below threshold
            % then sets all values outside this radius to 0.
            threshold = exp(-10);   % threshold amplitude for cutting off the greens function
            g_x = ifft2(obj.g0_k);
            g_radius = find(abs(g_x(1,:)) < threshold, 1);
            % checking grid size if green's function radius will fit
            % (otherwise new simulation is started with bigger sample)
            if (g_radius > obj.grid.y_margin) || (g_radius > obj.grid.x_margin)
                obj = wavesim(padarray(refractive_index, [g_radius, g_radius], n_center),bandwidth);
                return
            end            
            f_mask = @(x,y) x.^2+y.^2 < obj.grid.x_range(g_radius)^2;
            %%%g_x = g_x .* bsxfun(f_mask, obj.grid.x_range, obj.grid.y_range);
            %%%obj.g0_k = fft2(g_x);
            
            %% Potential map (V==k_r^2-k^2).
            obj.background = true(obj.grid.Ny,obj.grid.Nx);
            obj.background(obj.grid.y_margin+1:end-obj.grid.y_margin,obj.grid.x_margin+1:end-obj.grid.x_margin) = false;
            obj.V = ones(obj.grid.Ny, obj.grid.Nx) * (obj.k^2-obj.k_red^2);
            obj.V(~obj.background) = (refractive_index*2*pi/lambda).^2-obj.k_red^2;
            
            %% defining damping curve
            obj.g_curve = 1-linspace(0,1,g_radius).^2;
            damping_x = [ zeros(1,obj.grid.x_margin-g_radius), obj.g_curve(end:-1:1),...
                         ones(1,obj.grid.Nx - 2*obj.grid.x_margin),...
                         obj.g_curve, zeros(1,obj.grid.x_margin-g_radius)];
            
            damping_y = [ zeros(1,obj.grid.y_margin-g_radius), obj.g_curve(end:-1:1),...
                         ones(1,obj.grid.Ny - 2*obj.grid.y_margin),...
                         obj.g_curve, zeros(1,obj.grid.y_margin-g_radius)];        
        
            damping = damping_y' * damping_x;         
             
            %% Low pass filter potential function V
            obj.bandwidth = bandwidth;
            width = round(min(obj.grid.Nx, obj.grid.Ny)*bandwidth);
            w_x = [zeros(ceil((obj.grid.Nx-width)/2),1); tukeywin(width, 0.125); zeros(floor((obj.grid.Nx-width)/2),1)];
            w_y = [zeros(ceil((obj.grid.Ny-width)/2),1); tukeywin(width, 0.125); zeros(floor((obj.grid.Ny-width)/2),1)].';
            win2d = fftshift(bsxfun(@times, w_x, w_y));
            obj.V = ifft2(win2d.*fft2(obj.V));
                        
            obj.V = obj.V.*damping;
        end;
            
        function [E_x] = exec(obj, source)
            %% Execute simulation
            E_x = 0;
            energy_threshold = 1E-12;
            en=0;
            inter_step=5;
            source(obj.background) = 0;
                        
            tic;
            en_all = zeros(1,5000/inter_step);
            for it=1:5000
				%E_x = ifft2(obj.g0_k .* fft2(obj.V.*E_x+source));
                E_a = ifft2(obj.g0_k .* fft2(obj.V.*E_x+source));
				E_x = ifft2(conj(obj.g0_k) .* fft2(conj(obj.V).*(E_x-E_a))) + E_a;       
                if (mod(it,inter_step)==0)
                    toc;
                    en_prev = en;
                    en = wavesim.energy(E_x);
                    en_all(it/inter_step)=en;
                    plot(en_all); title(it); pause(0.1);
                    disp(['Added energy ', num2str(en-en_prev)]);
                    if (abs(en-en_prev) < energy_threshold)
                        disp(['Reached steady state in ' num2str(it) ' iterations']);
                        break;
                    end;
                    disp(['Ratio ', num2str(((en-en_prev)/en_prev)^(1/(2*inter_step)))]);
                    %plot(real(E_x(:,end/2)));title(it);pause(0.1);%imagesc(((real(E_x)))); colorbar; title(it); pause(0.1);
                    tic;
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
    end
end
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
            %% Constructing simgrid object and background mask            
            [rows,cols] = size(refractive_index);
            n_min = min(refractive_index(:));
            n_max = max(refractive_index(:));
            
            % padding of sample if uneven number of elements
            if mod(rows,2) == 1
                refractive_index(end+1,:) = sqrt(( n_max^2 + n_min^2 )/2);
                rows = rows + 1;
            end
            
            if mod(cols,2) == 1
                refractive_index(:,end+1) = sqrt(( n_max^2 + n_min^2 )/2);
                cols = cols + 1;
            end
                
            dx = 1/4;
            obj.grid = simgrid(cols,rows,dx);
            
            obj.background = true(obj.grid.Ny,obj.grid.Nx);
            obj.background(obj.grid.y_margin+1:end-obj.grid.y_margin,obj.grid.x_margin+1:end-obj.grid.x_margin) = false;

            %% Construct a wave simulation object for given refractive index map
            % currently, the refractive index map should have an average value of 1,
            % and the wavelength is fixed to 1 unit. 
            % window_size (step size) is the size of the exponential window
            % (in wavelengths)
            fftw('planner','patient'); %optimize fft2 and ifft2 at first use
            
            lambda = 1/dx;
            
                       
            obj.k = sqrt(( n_max^2 + n_min^2 )/2) * (2*pi/lambda);
            
            obj.epsilon = ( 2*pi/lambda )^2 * (n_max^2 - (n_max^2 + n_min^2 )/2) + 0.2;
            obj.k_red = sqrt(obj.k^2 + 1.0i*obj.epsilon); %obj.k+1.0i*obj.epsilon;
            
            disp(['epsilon = ',num2str(obj.epsilon)]);
      
            %% truncate Green function in x-space to have a finite support
            % for the convolution.
            % This operation is slighly
            % lossy, which means that the simulation itself is slightly lossy.
            % Usually, this slight loss is sufficient to make the divergence exponent
            % for the homogeneous case < 1
            
            % First find the distance at which g_x drops below threshold
            % then sets all values outside this radius to 0.
            f_g0_k = @(px, py) -1./(obj.k_red^2-px.^2-py.^2);
            obj.g0_k = bsxfun(f_g0_k, obj.grid.px_range, obj.grid.py_range);
            
            
            g_x = ifft2(obj.g0_k);
            obj.info.full_g0_k_max = max(abs(obj.g0_k(:)));
            obj.info.full_P = wavesim.energy(g_x);
                        
            % determining the green's function cut-off window
            threshold = exp(-10);   % threshold amplitude for cutting off the greens function
            g_radius = find(abs(g_x(1,:)) < threshold, 1);
            
            % checking grid size if green's function radius will fit
            % (otherwise new simulation is started with bigger sample)
            if g_radius > obj.grid.y_margin || g_radius > obj.grid.x_margin
                refractive_index_new = sqrt(( n_max^2 + n_min^2 )/2)*ones(rows + 2*g_radius, cols + 2*g_radius);
                refractive_index_new(g_radius+1:end-g_radius,g_radius+1:end-g_radius) = refractive_index;
                obj = wavesim(refractive_index_new,bandwidth);
                return
            end
            
            f_mask = @(x,y) x.^2+y.^2 < obj.grid.x_range(g_radius)^2;
            g_x = g_x .* bsxfun(f_mask, obj.grid.x_range, obj.grid.y_range);
            obj.g0_k = fft2(g_x);
            
            % Calculating the loss energy of the truncated green's function
            obj.info.trunc_g0_k_max = max(abs(obj.g0_k(:)));
            obj.info.trunc_P = wavesim.energy(g_x);
            obj.info.truncation_loss_P = obj.info.trunc_P/obj.info.full_P;
            obj.info.truncation_loss_g0_k_max = obj.info.trunc_g0_k_max/obj.info.full_g0_k_max;
            
            
            
            %% Potential map
            obj.V = ones(obj.grid.Ny, obj.grid.Nx);
            obj.V(obj.background) = - 1.0i*obj.epsilon;
            obj.V(~obj.background) = (refractive_index*2*pi/lambda).^2 - obj.k^2 - 1.0i*obj.epsilon;
            
            %% defining damping curve
%             % green function tail             
%             obj.g_curve = abs(g_x(1,1:g_radius));
%             obj.g_curve = obj.g_curve/max(obj.g_curve);
%             %obj.g_curve = 1 - obj.g_curve(end:-1:1)/max(obj.g_curve);
            
            obj.g_curve = 1-linspace(0,1,g_radius).^2;
            
            figure(5); plot(obj.g_curve);
            %% damping of the potential map (Method I: rectangular domain)
%             %damping of the potential map
            
            damping_x = [ zeros(1,obj.grid.x_margin-g_radius), obj.g_curve(end:-1:1),...
                         ones(1,obj.grid.Nx - 2*obj.grid.x_margin),...
                         obj.g_curve, zeros(1,obj.grid.x_margin-g_radius)];
            
            damping_y = [ zeros(1,obj.grid.y_margin-g_radius), obj.g_curve(end:-1:1),...
                         ones(1,obj.grid.Ny - 2*obj.grid.y_margin),...
                         obj.g_curve, zeros(1,obj.grid.y_margin-g_radius)];        
        
            damping = damping_y' * damping_x;         
            damping(damping < obj.g_curve(end-1)) = 0;
             
            %% Low pass filter potential function V
            
            obj.bandwidth = bandwidth;
            width = round(min(obj.grid.Nx, obj.grid.Ny)*bandwidth);
            w_x = [zeros(ceil((obj.grid.Nx-width)/2),1); tukeywin(width, 0.125); zeros(floor((obj.grid.Nx-width)/2),1)];
            w_y = [zeros(ceil((obj.grid.Ny-width)/2),1); tukeywin(width, 0.125); zeros(floor((obj.grid.Ny-width)/2),1)].';
            win2d = fftshift(bsxfun(@times, w_x, w_y));
            obj.V = ifft2(win2d.*fft2(obj.V));
                        
            %% applying damping to potential field
            figure(5); imagesc(real(damping)); colorbar;
            
            obj.V = obj.V.*damping;
        end;
            
        function [E_x] = exec(obj, source)
            %% Execute simulation
           
            %single pulse source
            E_x = full(source);
            energy_threshold = 1E-11;
            en=0;
            inter_step=10;
            
            source = zeros(obj.grid.Ny,obj.grid.Nx);
            
%             % cw-source
%             E_x = 0;
%             energy_threshold = 1E-12;
%             en=0;
%             inter_step=10;
%             source(obj.background) = 0;
            
            % predetermining green function convolution (new method II)
%             GS_k = obj.g0_k .*fft2( source );
%             VS_k = fft2(source.*conj(obj.V));            
            
            tic;
            for it=1:100
            %old method   
            %E_x = ifft2(obj.g0_k.*fft2(source-E_x.*obj.V));
                
                
                % New method I 
                E_x = ifft2(obj.g0_k .* (fft2(E_x.*obj.V+source) ...
                    - conj(obj.g0_k) .* fft2(E_x .* abs(obj.V).^2 + source.*conj(obj.V))));
                    + conj(obj.g0_k .* fft2(E_x.*conj(obj.V)));
                
                  % New method II
%                 B = fft2(E_x .* abs(obj.V).^2);
%                 C = fft2(E_x.*real(obj.V));
%                 D = 1i.*fft2(E_x.*imag(obj.V));
%                 
%                 E_x = ifft2(obj.g0_k .* ((C + D) ...
%                     - conj(obj.g0_k).* (B - VS_k)) ...
%                     + GS_k + conj(obj.g0_k) .* (C-D));
                
                if (mod(it,inter_step)==0)
                    toc;
                    en_prev = en;
                    en = wavesim.energy(E_x);
                    disp(['Added energy ', num2str(en-en_prev)]);
%                     if (abs(en-en_prev) < energy_threshold)
%                         disp(['Reached steady state in ' num2str(it) ' iterations']);
%                         break;
%                     end;
%                     disp(['Ratio ', num2str(((en-en_prev)/en_prev)^(1/(2*inter_step)))]);
                    figure(1); imagesc((abs(E_x))); colorbar; title(it); 
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
classdef wavesim
    %Simulation of the 2-D wave equation using a Born series approach
    % Ivo M. Vellekoop 2014
    
    properties
        V   % potential array used in simulation
        info % diagnostics information on convergence
        grid % simgrid object
     %%internal
        k % wave number
        k_red % k-i*epsilon
        epsilon % 1/'window size'
        g0_k% bare Green's function used in simulation
    end
        
    methods
        function obj=wavesim(grid, refractive_index, window_size, bandwidth)
            %% Construct a wave simulation object for given refractive index map
            % currently, the refractive index map should have an average value of 1,
            % and the wavelength is fixed to 1 unit. 
            % window_size (step size) is the size of the exponential window
            % (in wavelengths)
            lambda = 1;
            obj.grid = grid;
            obj.epsilon = 1/window_size;
            obj.k = 2*pi/lambda;
            obj.k_red = obj.k-1.0i*obj.epsilon;
            V_offset = obj.k_red^2-obj.k^2; %potential needed to convert from k-ieps^2 to k^2
            %% construct bare Green function.
            % The function is the solution for a medium with some
            % attenuation (k -> k-i epsilon).
            % todo: maybe an extra filter is needed to avoid aliasing
            % when the convolution of g0_k with fft2(V) wraps around p_max?
            % just low-pass filtering V may not be sufficient?
            %
            f_g0_k = @(px, py) 1./(obj.k_red^2-px.^2-py.^2);
            obj.g0_k = bsxfun(f_g0_k, grid.px_range, grid.py_range);
            
            %% truncate Green function in x-space to have a finite support
            % for the convolution.
            % This operation is slighly
            % lossy, which means that the simulation itself is slightly lossy.
            % Usually, this slight loss is sufficient to make the divergence exponent
            % for the homogeneous case < 1
            
            % First find the distance at which g_x drops below threshold
            % then sets all values outside this radius to 0.
            g_x = ifft2(obj.g0_k);
            obj.info.full_g0_k_max = max(abs(obj.g0_k(:)));
            obj.info.full_P = wavesim.energy(g_x);
            threshold = exp(-11);
            radius_index = find(abs(g_x(1,:))<threshold, 1);
            radius2 = grid.x_range(radius_index)^2;
            f_mask = @(x,y) x.^2+y.^2<radius2;
            g_x = g_x .* bsxfun(f_mask, grid.x_range, grid.y_range);
            obj.g0_k = fft2(g_x);
            obj.info.trunc_g0_k_max = max(abs(obj.g0_k(:)));
            obj.info.trunc_P = wavesim.energy(g_x);
            obj.info.truncation_loss_P = obj.info.trunc_P/obj.info.full_P;
            obj.info.truncation_loss_g0_k_max = obj.info.full_g0_k_max/obj.info.full_g0_k_max;
            
            %% Set up refractive index map for simulation. 
            % The map should have an average refractive index of 1 (or close?)
            % The map is automatically scaled to size Nx-2B, Ny-2B, where
            % B is the boundar y width. The boundary width defaults to 1.5
            % times the window size.
            B = radius_index*2;
            obj.V = ones(grid.Ny, grid.Nx); %background refractive index
            obj.V(B+1:end-B, B+1:end-B) = imresize(refractive_index, [grid.Ny-2*B, grid.Nx-2*B]);
            obj.V = (obj.V.^2-1) * obj.k^2 + V_offset; % convert refractive index to potential (remove '0-potential', the potential for which g0_k was constructed) 
        
            %% Construct boundaries. The potential in the boundary slowly
            % varies from lossless propagation (k->V=k_red^2-k^2) to damped propagation
            % (k_red->V=0)
            %
            % use quadratic increase from epsilon = 0
            epsilon_boundary = V_offset*[1-linspace(0,1,B/2).^2, zeros(1, B/2)];%(obj.k-1.0i*(1:B).^2/B^2*obj.epsilon).^2 - obj.k^2 + V_offset;
%            epsilon_x = [fliplr(epsilon_boundary), zeros(1, grid.Nx-2*B), epsilon_boundary];
 %           epsilon_y = [fliplr(epsilon_boundary), zeros(1, grid.Ny-2*B), epsilon_boundary].';
  %          epsilon_xy = (obj.k-1.0i*@bsxfun(@plus, epsilon_x, epsilon_y)).^2-obj.k_red^2;
            obj.V(1:B, :) = fliplr(epsilon_boundary).' * ones(1, grid.Nx);
            obj.V(end-B+1:end, :) = epsilon_boundary.' * ones(1, grid.Nx);
            obj.V(B+1:end-B, 1:B) = 100;
            obj.V(B+1:end-B, end-B+1:end) = 100;
            obj.V(:, 1:B) = min(obj.V(:, 1:B), ones(grid.Ny, 1) * fliplr(epsilon_boundary));
            obj.V(:, end-B+1:end) = min(obj.V(:,end-B+1:end), ones(grid.Ny, 1) * epsilon_boundary);
  %        keyboard;
            %% Low pass filter potential function V
            width = round(min(grid.Nx, grid.Ny)*bandwidth);
            w_x = [zeros(ceil((grid.Nx-width)/2),1); tukeywin(width, 0.125); zeros(floor((grid.Nx-width)/2),1)];
            w_y = [zeros(ceil((grid.Ny-width)/2),1); tukeywin(width, 0.125); zeros(floor((grid.Ny-width)/2),1)].';
            win2d = fftshift(bsxfun(@times, w_x, w_y));
            obj.V = ifft2(win2d.*fft2(obj.V));
        end;
        function E_tot = exec(obj, source)
            %% Execute simulation
            E_tot = 0;
            E_x = 0;
            threshold = exp(-35);
            energy_threshold = 1E-9;
            en_prev=1;
            inter_step=5;
            for it=1:1000
                if it<2
                    E_x = ifft2(obj.g0_k.*fft2(E_x.*obj.V+source));
                else
                    E_x = ifft2(obj.g0_k.*fft2(E_x.*obj.V));
                end;
                
                %apply threshold to avoid accumulation of machine precision errors
                E_x = E_x .* (abs(E_x)>threshold);
                E_tot = E_tot + E_x;

                if (mod(it,inter_step)==1)
                    %dispField = bsxfun(@multiply, E_x, reference'); %interesting image: shows scattering
                    %mean free path?
                    %imagesc(real(E_x));
                    %imagesc(log(abs(fft2(E_x)))); colorbar;
                    imagesc(log(abs(E_x))); colorbar;
                    %plot(log(abs(E_x(:,end/2))));
                    %plot(p_range/k, abs(fft(E_x(:,end/2))).^2);
                    %rgb = imoverlay(mat2gray(abs(scattered)), t<0, [0 1 0], 0.25);
                    %imshow(rgb, 'InitialMagnification', 'fit');
                    %en = energy(E_x);
                    en = 2*real(E_tot(:)'*E_x(:));%energy increase for this step
                    disp(['Added energy ', num2str(en)]);
                    if (abs(en) < energy_threshold)
                        disp('Reached steady state');
                        break;
                    end;
                    disp(['Ratio ', num2str((en/en_prev)^(1/(2*inter_step)))]);
                    en_prev = en;
                    title(it); pause(0.5);
                end;
            end;
            imagesc(log(abs(E_tot))); colorbar;
        end;
    end
    methods(Static)
        function en = energy(E_x)
            en= E_x(:)' * E_x(:);
        end;
    end
end


classdef Medium
    % Medium - Generates a sample object for use in wave
    % simulations (wavesim, PSTD)
    % Ivo Vellekoop & Gerwin Osnabrugge 2016-2020
    properties
        e_r             % relative dielectric constand map with boundaries appended
        e_r_max         % maximum real part of obj.e_r
        e_r_min         % minimum real part of obj.e_r
        e_r_center      % (e_r_max+e_r_min)/2
        boundary_type;  % type of boundary layer used around sample
        Bl;             % left boundary widths (absorbing boundaries + padding)
        Br;             % right boundary widths (absorbing boundaries + padding)
        n_media;        % number of potential maps stored in memory (multiple 
                        % submedia are used for anti-aliasing algorithm)
    end
    methods
        function [obj,grid] = Medium(refractive_index, options)
            % Medium - Generates a sample object for use in wave
            % simulations (wavesim, PSTD)
            %
            % internally, the object stores a map of the relative dielectric
            % constant (e_r), which is the refractive index squared. The e_r map
            % is padded with boundaries, and then expanded to the next
            % multiple of 2^N in each dimension (for fast fourier transform)
            %
            % refractive_index   = refractive index map, may be complex, need not
            %                      be square. Can be 2-D or 3-D.
            % options.pixel_size = size of a grid pixel in any desired unit (e.g.
            % microns)
            % options.boundary_widths = vector with widths of the absorbing
            %                          boundary (in pixels) for each dimension. Set element
            %                          to 0 for periodic boundary.
            % options.boundary_strength = maximum value of imaginary part of e_r
            %                             in the boundary (for 'PBL' only)
            % options.boundary_type = boundary type. Currently supports 'window' (default) and
            % 'PBL1-5'
            % options.ar_width = width of anti-reflection layer (for
            % 'window' only). Should be < boundary_widths. When omitted, a
            % default fraction of the boundary width will be used
            % (pre-factor may change in the future)
            %
            % New feature:
            % options.medium_wiggle (3x1 boolean array, corresponding with
            % dimensions y,x,z).
            % when enabled the medium will be subdivided into multiple sub
            % media for anti-aliasing
            
            %% check validity of input sample
            if min(imag(refractive_index(:))) < 0
                error('Medium cannot have gain, imaginary part of refractive index should be negative');
            end
            
            %% Set default values and check validity of inputs
            assert(numel(refractive_index) >= 1);
            assert(numel(options.boundary_widths) >= ndims(refractive_index));
            assert(numel(options.boundary_widths) <= 3);
            if ~isfield(options, 'boundary_type')
                options.boundary_type = 'window';
                obj.boundary_type = options.boundary_type;
            else    
                obj.boundary_type = 'window'; %new default boundary type!
            end       
            if ~isfield(options, 'medium_wiggle')
                medium_wiggle = false; % by default medium wiggle is disabled in all dimensions
            end
            
            %% calculate e_r and min/max values
            e_r = refractive_index.^2;
            obj.e_r_min = min(real(e_r(:)));
            obj.e_r_max = max(real(e_r(:)));
            obj.e_r_center = (obj.e_r_min + obj.e_r_max)/2;
            
            % subsample medium into smaller media when medium_wiggle is
            % enabled
            [obj.e_r,obj.n_media] = Medium.subsample(e_r,medium_wiggle);
            
            %% construct coordinate set and calculate padding width
            % padds to next efficient size for fft in each dimension, and makes sure to
            % append at least 'boundary_widths' pixels on both sides.
            sz = Simulation.make3(size(obj.e_r{1}), 1);
            bw = Simulation.make3(options.boundary_widths * 2, 0);
            grid = SimulationGrid(sz + bw, options.pixel_size, bw==0);
            
            % calculate effective boundary width (absorbing boundaries + padding)
            % on left and right hand side, respectively.
            obj.Bl = ceil((grid.N - sz) / 2);
            obj.Br = floor((grid.N - sz) / 2);
            
            %% add boundary layer around the medium
            for i_medium = 1:obj.n_media                
                obj.e_r{i_medium} = obj.add_boundary_layer(obj.e_r{i_medium}, options);
            end
        end
    end
    methods(Access=protected)     
        function e_r = add_boundary_layer(obj, e_r, options)        
            % adds a boundary layer around the medium
            e_r = Medium.extrapolate(e_r, obj.Bl, obj.Br);
            if (strcmp(options.boundary_type(1:3), 'PBL') && any(obj.Bl~=0)) 
                % add absorption to the boundary layer (deprecated) 
                e_r = obj.polynomial_boundary_layer(e_r, options);
            end
        end
        function e_r = polynomial_boundary_layer(obj, e_r, options)
            % the polynomial boundary layer adds absorption in such a way 
            % to minimize  reflection of a normally incident beam
            warning('PBL boundary conditions are deprecated and may be removed in the future: use window');
            
            %the shape of the boundary is determined by f_boundary_curve, a function
            %that takes a position (in pixels, 0=start of boundary) and returns
            %Delta e_r for the boundary.
            %todo: e_0 per row/column? or per side?
            e_0 = mean(e_r(:));
            k0 = sqrt(e_0)*2*pi/ (options.lambda / options.pixel_size); %k0 in 1/pixels
            % maximum value of the boundary (see Mathematica file = c(c-2ik0) = boundary_strength)
            % ||^2 = |c|^2 (|c|^2 + 4k0^2)   [assuming c=real, better possible when c complex?]
            % when |c| << |2k0| we can approximage: boundary_strength = 2 k0 c
            c = options.boundary_strength*k0^2 / (2*k0);
            switch (options.boundary_type)
                %case 'PBL' %Nth order smooth?
                %    N=3;
                %    f_boundary_curve = @(r) 1/k0^2*(c^(N+2)*r.^N.*(N+1.0+(2.0i*k0-c)*r)) ./ (factorial(N)*exp(c*r));
                case 'PBL5' %5th order smooth
                    f_boundary_curve = @(r) 1/k0^2*(c^7*r.^5.*(6.0+(2.0i*k0-c)*r)) ./ (720+720*c*r+360*c^2*r.^2+120*c^3*r.^3+30*c^4*r.^4+6*c^5*r.^5+c^6*r.^6);
                case 'PBL4' %4th order smooth
                    f_boundary_curve = @(r) 1/k0^2*(c^6*r.^4.*(5.0+(2.0i*k0-c)*r)) ./ (120+120*c*r+60*c^2*r.^2+20*c^3*r.^3+5*c^4*r.^4+c^5*r.^5);
                case 'PBL3' %3rd order smooth
                    f_boundary_curve = @(r) 1/k0^2*(c^5*r.^3.*(4.0+(2.0i*k0-c)*r)) ./ (24+24*c*r+12*c^2*r.^2+4*c^3*r.^3+c^4*r.^4);
                case 'PBL2' %2nd order smooth
                    f_boundary_curve = @(r) 1/k0^2*(c^4*r.^2.*(3.0+(2.0i*k0-c)*r)) ./ (6+6*c*r+3*c^2*r.^2+c^3*r.^3);
                case 'PBL1' %1st order smooth
                    f_boundary_curve = @(r) 1/k0^2*(c^3*r.*(2.0+(2.0i*k0-c)*r)) ./ (2.0+2.0*c*r+c^2*r.^2) / k0^2; %(divide by k0^2 to get relative e_r)
                otherwise
                    error(['unknown boundary type' options.boundary_type]);
            end
            roi_size = Simulation.make3(size(e_r), 1) - obj.Bl - obj.Br;
            x = [(obj.Bl(2):-1:1), zeros(1, roi_size(2)), (1:obj.Br(2))];
            y = [(obj.Bl(1):-1:1), zeros(1, roi_size(1)), (1:obj.Br(1))];
            z = [(obj.Bl(3):-1:1), zeros(1, roi_size(3)), (1:obj.Br(3))];
            x = reshape(x, [1,length(x),1]);
            y = reshape(y, [length(y),1,1]);
            z = reshape(z, [1,1,length(z)]);
            e_r = e_r + f_boundary_curve(sqrt(x.^2 + y.^2 + z.^2));
        end   
    end
    methods (Static)
        function e_r = extrapolate(e_r, Bl, Br)
            % pads the permittivity map e_r to total boundary widths Bl and
            % Br. the new pixels will be filled with repeated edge pixels
            % the boundaries are added to both sides.
            e_r = padarray(e_r, Bl, 'replicate', 'pre');
            e_r = padarray(e_r, Br, 'replicate', 'post');
        end
        function [e_r_set,n_media] = subsample(e_r, medium_wiggle)
            % function used to subsample medium into smaller submedia for
            % anti-aliasing. Odd indices represent grid points on a
            % right-shifted grid and even indices represent grid points on a
            % left-shift grid.
            
            Nx = [size(e_r,1), size(e_r,2), size(e_r,3)]; % size sample
            
            % check medium_wiggle input
            if medium_wiggle == true
                medium_wiggle = Nx > 1;
            elseif medium_wiggle == false
                medium_wiggle = false(1,3);
            elseif numel(medium_wiggle) < 3
                medium_wiggle(end+1:3) = false;
            end
       
            % check if all wiggled direction have an even number of 
            % gridpoints. If not, append to make grid even
            append = double(mod(Nx,2) & medium_wiggle);
            e_r = padarray(e_r, append, 'replicate', 'post');

            % determine wiggle directions [y;x;z]
            medium_wiggles = wiggle_perm(medium_wiggle);
            
            % generate wiggled submedia 
            n_media = 2^(sum(medium_wiggle));
            e_r_set = cell(n_media,1);
            for i_medium = 1:n_media
                wshift = (medium_wiggles(:,i_medium) == 1);
                e_r_set{i_medium} = e_r(1+wshift(1):1+medium_wiggle(1):end,...
                    1+wshift(2):1+medium_wiggle(2):end,...
                    1+wshift(3):1+medium_wiggle(3):end);
            end
        end
        
        %%% testing methods
        function test_subsample()
            % unit test for Medium.subsample method
            A = 1+rand(10,10,10);
            
            % test 1
            disp('test 1: no medium wiggle enabled...');
            [B,n_media] = Medium.subsample(A,[0,0,0]);
            B_expected{1} = A;
            if isequal(B{1},B_expected{1}) && n_media == 1
                disp('passed');
            else
                disp('failed');
            end
            
            % test 2
            disp('test 2: single medium wiggle...');
            [B,n_media] = Medium.subsample(A,[1,0,0]);
            B_expected{1} = A(2:2:end,:,:);
            B_expected{2} = A(1:2:end,:,:);
            
            for i_medium = 1:numel(B)
                if ~isequal(B{i_medium},B_expected{i_medium})
                    disp('failed');
                    break;
                end
                if i_medium == numel(B) && n_media == 2
                    disp('passed');
                end
            end
            
            % test 3:
            disp('test 3: medium wiggle in 2D...');
            [B,n_media] = Medium.subsample(A,[1,0,1]);
            B_expected{1} = A(2:2:end,:,2:2:end);
            B_expected{2} = A(1:2:end,:,2:2:end);
            B_expected{3} = A(2:2:end,:,1:2:end);
            B_expected{4} = A(1:2:end,:,1:2:end);

            for i_medium = 1:numel(B)
                if ~isequal(B{i_medium},B_expected{i_medium}) 
                    disp('failed');
                    break;
                end
                if i_medium == numel(B) && n_media == 4
                    disp('passed');
                end
            end
            
            % test 4:
            disp('test 4: medium wiggle in 3D...');
            [B,n_media] = Medium.subsample(A,[1,1,1]);
            B_expected{1} = A(2:2:end,2:2:end,2:2:end);
            B_expected{2} = A(1:2:end,2:2:end,2:2:end);
            B_expected{3} = A(2:2:end,1:2:end,2:2:end);
            B_expected{4} = A(1:2:end,1:2:end,2:2:end);
            B_expected{5} = A(2:2:end,2:2:end,1:2:end);
            B_expected{6} = A(1:2:end,2:2:end,1:2:end);
            B_expected{7} = A(2:2:end,1:2:end,1:2:end);
            B_expected{8} = A(1:2:end,1:2:end,1:2:end);
            for i_medium = 1:numel(B)
                if ~isequal(B{i_medium},B_expected{i_medium}) 
                    disp('failed');
                    break;
                end
                if i_medium == numel(B) && n_media == 8
                    disp('passed');
                end
            end
            
            % test 5:
            disp('test 5: test sample with odd number of grid points...');
            A = rand(11,11,1);
            [B,~] = Medium.subsample(A,[0,1]);
            B_expected{1} = [A(:,2:2:end),A(:,end)];   % appended direction
            B_expected{2} = A(:,1:2:end);

            if isequal(B{1},B_expected{1}) && isequal(B{2},B_expected{2})
                disp('passed');
            else
                disp('failed');
            end
        end
    end
end
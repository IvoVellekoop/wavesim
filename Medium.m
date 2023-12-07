classdef Medium
    % Medium - Generates a sample object for use in wave
    % simulations (wavesim, PSTD)
    % Ivo Vellekoop & Gerwin Osnabrugge 2016-2020
    properties
        boundary_type = 'ARL';  % type of boundary layer used around sample (default: anti-reflection layer)
        e_r             % relative dielectric constand map with added boundary layers
        e_r_max         % maximum real part of obj.e_r
        e_r_min         % minimum real part of obj.e_r
        e_r_center      % (e_r_max+e_r_min)/2
        Bl;             % left boundary widths (absorbing boundaries + padding)
        Br;             % right boundary widths (absorbing boundaries + padding)
    end
    methods
        function [obj,grid] = Medium(refractive_index, options)
            % Medium - Generates a sample object for use in wave
            % simulations (wavesim, PSTD)
            %
            % internally, the object stores a map of the relative dielectric
            % constant (e_r), which is the refractive index squared. The 
            % e_r map is padded with boundaries, and then expanded in each 
            % dimension to an efficient size for fast fourier transform 
            %
            % refractive_index   = refractive index map, may be complex, 
            %                      need not be square. Can be 2-D or 3-D.
            % options.pixel_size = size of a grid pixel in any desired unit .
            %                      (e.g microns)
            % options.boundary_widths = vector with widths of the absorbing
            %                           boundary (in pixels) for each dimension. 
            %                           Set element to 0 for periodic boundary.
            % options.boundary_strength = maximum value of imaginary part 
            %                             of e_r in the boundary (for 'PBL' only)
            % options.boundary_type = boundary type. Currently supports 
            %                         'ARL' (default) and 'PBL1-5'
            
            %% check validity of input sample
            if min(imag(refractive_index(:))) < 0
                error('Medium cannot have gain, imaginary part of refractive index should be positive');
            end
            
            % set default boundary type
            if isfield(options, 'boundary_type')
                obj.boundary_type = options.boundary_type;
            end
            
            %% Set default values and check validity of inputs
            assert(numel(refractive_index) >= 1); % not sure why this is needed
            assert(numel(options.boundary_widths) >= ndims(refractive_index), "boundary_widths array is too short");
            assert(numel(options.boundary_widths) <= 3, "boundary_widths array is too long");
            assert(all(round(options.boundary_widths)==options.boundary_widths), "boundary_widths should be integer numbers")
            
            %% calculate e_r and min/max values
            e_r = refractive_index.^2;
            obj.e_r_min = min(real(e_r(:)));
            obj.e_r_max = max(real(e_r(:)));
            obj.e_r_center = (obj.e_r_min + obj.e_r_max)/2;
            
            %% construct coordinate set and calculate padding width
            % padds to next efficient size for fft in each dimension, and makes sure to
            % append at least 'boundary_widths' pixels on both sides.
            sz = Simulation.make3(size(e_r), 1);
            bw = Simulation.make3(options.boundary_widths * 2, 0);
            grid = SimulationGrid(sz + bw, options.pixel_size, bw==0);
            
            % calculate effective boundary width (absorbing boundaries + padding)
            % on left and right hand side, respectively.
            obj.Bl = ceil((grid.N - sz) / 2);
            obj.Br = floor((grid.N - sz) / 2);
            
            %% add boundary layer around the medium             
            obj.e_r = obj.add_boundary_layer(e_r, options);
        end
    end
    methods(Access=protected)     
        function e_r = add_boundary_layer(obj, e_r, options)        
            % adds a boundary layer around the medium
            e_r = Medium.extrapolate(e_r, obj.Bl, obj.Br);
            if (strcmp(obj.boundary_type(1:3), 'PBL') && any(obj.Bl~=0)) 
                % add absorption to the boundary layer (deprecated) 
                e_r = obj.polynomial_boundary_layer(e_r, options);
            end
        end
        function e_r = polynomial_boundary_layer(obj, e_r, options)
            % the polynomial boundary layer adds absorption in such a way 
            % to minimize  reflection of a normally incident beam
            warning('PBL boundary conditions are deprecated and may be removed in the future: use ARL');
            
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
    end
end
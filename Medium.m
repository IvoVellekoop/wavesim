classdef Medium
% Medium - Generates a sample object for use in wave
% simulations (wavesim, PSTD)
% Ivo Vellekoop 2017
    properties
        e_r % relative dielectric constand map with boundaries appended
        e_r_max % maximum real part of obj.e_r
        e_r_min % minimum real part of obj.e_r
        e_r_center % (e_r_max+e_r_min)/2
        roi % size of the original refractive index map without padding
        grid %  simgrid object, with x and k ranges
        filters = []% profiles for windowed absorbing boundaries
        leakage = 0
    end
    methods
        function obj = Medium(refractive_index, options)
            % Medium - Generates a sample object for use in wave
            % simulations (wavesim, PSTD, FDTD)
            % 
            % internally, the object stores a map of the relative dielectric
            % constant (e_r), which is the refractive index squared. The e_r map
            % is padded with absorbing boundaries, and then expanded to the next
            % multiple of 2^N in each dimension (for fast fourier transform)
            %
            % refractive_index   = refractive index map, may be complex, need not
            %                      be square. Can be 2-D or 3-D.
            % options.pixel_size = size of a grid pixel in any desired unit (e.g.
            % microns)
            % options.lambda     = wavelength (in same units as pixel_size)
            % options.boundary_widths = vector with widths of the absorbing
            %                          boundary (in pixels) for each dimension. Set element
            %                          to 0 for periodic boundary.
            % options.boundary_strength = maximum value of imaginary part of e_r
            %                             in the boundary (for 'PML' only)
            % options.boundary_type = boundary type. Currently supports 'window' (default) and
            % 'PML1-5'
            % options.ar_width = width of anti-reflection layer (for
            % 'window' only). Should be < boundary_widths. When omitted, a
            % default fraction of the boundary width will be used
            % (pre-factor may change in the future)
            %
            %% Set default values and check validity of inputs
            assert(numel(options.boundary_widths) == ndims(refractive_index));
            if ~isfield(options, 'boundary_type')
                options.boundary_type = 'window'; %new default boundary type!
            end
            if ~isfield(options, 'ar_width')
                options.ar_width = round(options.boundary_widths*0.75); %this factor may change into a more optimal one!!!
            end
            
            %% calculate e_r and min/max values
            obj.e_r = refractive_index.^2;
            obj.e_r_min = min(real(obj.e_r(:)));
            obj.e_r_max = max(real(obj.e_r(:)));
            obj.e_r_center = (obj.e_r_min + obj.e_r_max)/2;
            
            % construct coordinate set. 
            % padds to next efficient size for fft in each dimension, and makes sure to
            % append at least 'boundary_widths' pixels on both sides.
            sz = Medium.make3(size(obj.e_r), 1);
            bw = Medium.make3(options.boundary_widths * 2, 0);
            obj.grid = SimulationGrid(sz + bw, options.pixel_size, bw==0);

            % applies the padding, extrapolating the refractive index map to
            % into the added regions
            [obj.e_r, obj.roi, Bl, Br] = Medium.extrapolate(obj.e_r, obj.grid.N);
            [obj.e_r, obj.leakage] = Medium.add_absorbing_boundaries(obj.e_r, Bl, Br, options); 
            obj.filters = Medium.edge_filters(obj.grid.N, Bl, Br, options);
        end
    end
    methods (Static)
        function sz = make3(sz, pad)
            % makes sure sz is a row vector with exactly 3 elements
            % (pads when needed)
            % (this function is needed because 'size' by default removes
            %  trailing singleton dimensions)
            sz = sz(:).'; %make row vector
            if numel(sz) < 3
                sz((end+1):3) = pad;
            end
        end
        
        function [e_r, roi, Bl, Br] = extrapolate(e_r, new_size)
            %% Expands the permittivity map to 'new_size'
            % The new pixels will be filled with repeated edge pixels
            
            % Calculate effective boundary width (absorbing boundaries + padding)
            % on left and right hand side, respectively.
            roi_size = Medium.make3(size(e_r), 1);
            Bl = ceil((new_size - roi_size) / 2); 
            Br = floor((new_size - roi_size) / 2); %effective boundary width (absorbing boundaries + padding)
            
            % the boundaries are added to both sides. Remember where is the region of interest
            roi = [Bl + 1; Bl+roi_size];
            e_r = padarray(e_r, Bl, 'replicate', 'pre');
            e_r = padarray(e_r, Br, 'replicate', 'post');
        end
        
        function [e_r, leakage] = add_absorbing_boundaries(e_r, Bl, Br, options)
            %only for (now deprecated) PML boundary conditions:
            %Adds absorption in such a way to minimize reflection of a 
            %normally incident beam
            %
            if (~strcmp(options.boundary_type(1:3), 'PML') || all(Bl==0)) 
                leakage = [];
                return;
            end
            warning('PML boundary conditions are deprecated and may be removed in the future: use window');
            %the shape of the boundary is determined by f_boundary_curve, a function
            %that takes a position (in pixels, 0=start of boundary) and returns
            %Delta e_r for the boundary. 
            Bmax = max(Br); %used to calculate expected amount of leakage through boundary
            %todo: e_0 per row/column? or per side?
            e_0 = mean(e_r(:));
            k0 = sqrt(e_0)*2*pi/ (options.lambda / options.pixel_size); %k0 in 1/pixels
            % maximum value of the boundary (see Mathematica file = c(c-2ik0) = boundary_strength)
            % ||^2 = |c|^2 (|c|^2 + 4k0^2)   [assuming c=real, better possible when c complex?]
            % when |c| << |2k0| we can approximage: boundary_strength = 2 k0 c
            c = options.boundary_strength*k0^2 / (2*k0);
            switch (options.boundary_type)
                %case 'PML' %Nth order smooth?
                %    N=3;
                %    f_boundary_curve = @(r) 1/k0^2*(c^(N+2)*r.^N.*(N+1.0+(2.0i*k0-c)*r)) ./ (factorial(N)*exp(c*r));
                %    obj.leakage = exp(-c*Bmax)*exp(c*Bmax);
                case 'PML5' %5th order smooth
                    f_boundary_curve = @(r) 1/k0^2*(c^7*r.^5.*(6.0+(2.0i*k0-c)*r)) ./ (720+720*c*r+360*c^2*r.^2+120*c^3*r.^3+30*c^4*r.^4+6*c^5*r.^5+c^6*r.^6);
                    leakage = exp(-c*Bmax)*(720+720*c*Bmax+360*c^2*Bmax.^2+120*c^3*Bmax.^3+30*c^4*Bmax.^4+6*c^5*Bmax.^5+c^6*Bmax.^6)/24;
                case 'PML4' %4th order smooth
                    f_boundary_curve = @(r) 1/k0^2*(c^6*r.^4.*(5.0+(2.0i*k0-c)*r)) ./ (120+120*c*r+60*c^2*r.^2+20*c^3*r.^3+5*c^4*r.^4+c^5*r.^5);
                    leakage = exp(-c*Bmax)*(120+120*c*Bmax+60*c^2*Bmax.^2+20*c^3*Bmax.^3+5*c^4*Bmax.^4+c^5*Bmax.^5)/24;
                case 'PML3' %3rd order smooth
                    f_boundary_curve = @(r) 1/k0^2*(c^5*r.^3.*(4.0+(2.0i*k0-c)*r)) ./ (24+24*c*r+12*c^2*r.^2+4*c^3*r.^3+c^4*r.^4);
                    leakage = exp(-c*Bmax)*(24+24*c*Bmax+12*c^2*Bmax.^2+4*c^3*Bmax.^3+c^4*Bmax.^4)/24;
                case 'PML2' %2nd order smooth
                    f_boundary_curve = @(r) 1/k0^2*(c^4*r.^2.*(3.0+(2.0i*k0-c)*r)) ./ (6+6*c*r+3*c^2*r.^2+c^3*r.^3);
                    leakage = exp(-c*Bmax)*(6+6*c*Bmax+3*c^2*Bmax.^2+c^3*Bmax.^3)/6;
                case 'PML1' %1st order smooth
                    f_boundary_curve = @(r) 1/k0^2*(c^3*r.*(2.0+(2.0i*k0-c)*r)) ./ (2.0+2.0*c*r+c^2*r.^2) / k0^2; %(divide by k0^2 to get relative e_r)
                    leakage = exp(-c*Bmax)*(2+2*c*Bmax+c^2*Bmax.^2)/2;
                otherwise
                    error(['unknown boundary type' obj.boundary_type]);
            end
            roi_size = Medium.make3(size(e_r), 1) - Bl - Br;
            x = [(Bl(2):-1:1), zeros(1, roi_size(2)), (1:Br(2))];
            y = [(Bl(1):-1:1), zeros(1, roi_size(1)), (1:Br(1))];
            z = [(Bl(3):-1:1), zeros(1, roi_size(3)), (1:Br(3))];
            x = reshape(x, [1,length(x),1]);
            y = reshape(y, [length(y),1,1]);
            z = reshape(z, [1,1,length(z)]);
            e_r = e_r + f_boundary_curve(sqrt(x.^2 + y.^2 + z.^2));
        end
        
        function filters = edge_filters(full_size, Bl, Br, options)
            % Note: currently experimental, may be removed.
            % construct filters that simply nullify the field outside the roi
            % (but utilize a smooth transition function)
            filters = cell(3,1);
            if (~strcmp(options.boundary_type, 'window') || all(Bl==0)) 
                return;
            end
            
            for dim=1:3
                bl = Bl(dim); %width of added boundary
                br = Br(dim);
                roi_size = full_size(dim) - bl - br;
                if bl > 0
                    L = options.ar_width(dim); % width of window
                   % window = parzenwin(L);
                    window = nuttallwin(L);
                    smoothstep = cumsum(window)/sum(window); % integrate window to get step function
                    filt = [zeros(bl-L, 1); smoothstep; ones(roi_size, 1); flipud(smoothstep); zeros(br-L, 1)];
                    filters{dim} = reshape(filt, circshift([1, 1, length(filt)], [0,dim]));
                end
            end
        end
    end
end
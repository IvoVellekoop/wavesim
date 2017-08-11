classdef simgrid
    %grid parameters for simulation
    % Ivo M. Vellekoop
    properties
        N  % number of grid points in y,x,z-dimension
        dx  % grid resolution (x and y)
        padding % amount of added zero padding
        x_range % grid point coordinates in x-dimension
        y_range % grid point coordinates in y-dimension
        z_range % grid point coordinates in z-dimension
        px_range % grid point coordinates in Fourier transformed x-dimension
        py_range % grid point coordinates in Fourier transformed y-dimension
        pz_range % grid point coordinates in Fourier transformed z-dimension
    end
    methods
       function obj = simgrid(min_size, dx)
            % Construct a wave simulation grid object with specified width,
            % min_size = [height,width], minimum required size (will be rounded up)
			% dx = step size of grid (in arbitrary units)
			%% setup coordinates
            if numel(min_size) ~= 3
                error ('Only 3-D structures are supported');
            end
            N = pow2(nextpow2(min_size));
            obj.padding = N - min_size; %total amoung of zero padding. Usually placed at right and bottom sides only (non-centric)
            obj.dx = dx;
            obj.x_range = (0:(N(2)-1))*dx; %obj.x_range = linspace(min_x,max_x,obj.N(2));
            obj.y_range = (0:(N(1)-1)).'*dx; %obj.y_range = linspace(min_y,max_y,obj.N(1));
            obj.z_range = (0:(N(3)-1))*dx; %obj.z_range = linspace(min_x,max_x,obj.N(2));
            obj.px_range = 2*pi*simgrid.symrange(N(2))/(dx*N(2));
            obj.py_range = 2*pi*simgrid.symrange(N(1)).'/(dx*N(1));
            obj.pz_range = reshape(2*pi*simgrid.symrange(N(3))/(dx*N(3)), [1, 1, N(3)]);
            obj.N = N;
       end
       function retval = p2(obj)
          %returns a 2 or 3 dimensional array with the magnitude of p squared.
          retval = simgrid.dist2_3d(obj.px_range, obj.py_range, obj.pz_range);
       end
    end
	methods(Static)
        function retval = dist2_2d(x, y)
            %returns a 2 dimensional array with x^2+y^2
            f_p2 = @(xi, yi) xi.^2 + yi.^2;
            retval = bsxfun(f_p2, x, y);
        end
        function retval = dist2_3d(x, y, z)
            %returns a 3 dimensional array with x^2+y^2+z^2
             f_pz = @(xy2, zi) xy2 + zi.^2;
             retval = bsxfun(f_pz, simgrid.dist2_2d(x, y), z);
        end
        function range = symrange(N)
            if N==1
                range = 0;
            else
                range = fftshift(-N/2:N/2-1);
            end
        end
    end
end


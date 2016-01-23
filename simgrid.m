classdef simgrid
    %grid parameters for simulation
    % Ivo M. Vellekoop
    properties
        N  % number of grid points in y,x,z-dimension
        Nred % same as N, but leaving out singleton dimensions
        dx  % grid resolution (x and y)
        padding % amount of added zero padding
        x_range % grid point coordinates in x-dimension
        y_range % grid point coordinates in y-dimension
        z_range % grid point coordinates in z-dimension
        px_range % grid point coordinates in Fourier transformed x-dimension
        py_range % grid point coordinates in Fourier transformed y-dimension
        pz_range % grid point coordinates in Fourier transformed z-dimension
        dimension % dimension of the simulation (= length of min_size in constructor)
    end
    methods
       function obj = simgrid(min_size, dx)
            % Construct a wave simulation grid object with specified width,
            % min_size = [height,width], minimum required size (will be rounded up)
			% dx = step size of grid (in arbitrary units)
			%% setup coordinates
            obj.dimension = length(min_size);
            if (obj.dimension == 1)
                min_size(2:3) = 1; % extend 1-D to 3-D
            elseif (obj.dimension == 2)
                min_size(3) = 1; %extend 2-D to 3-D
            elseif (obj.dimension ~= 3)
                error (['dimension ' num2str(obj.dimension) ' is not supported.']);
            end;
            obj.N = pow2(nextpow2(min_size));
            obj.Nred = obj.N(1:obj.dimension); %leave out singleton dimensions
            obj.padding = obj.N - min_size; %total amoung of zero padding. Usually placed at right and bottom sides only (non-centric)
            obj.dx = dx;
            obj.x_range = (0:(obj.N(2)-1))*dx; %obj.x_range = linspace(min_x,max_x,obj.N(2));
            obj.y_range = (0:(obj.N(1)-1)).'*dx; %obj.y_range = linspace(min_y,max_y,obj.N(1));
            obj.z_range = (0:(obj.N(3)-1))*dx; %obj.z_range = linspace(min_x,max_x,obj.N(2));
            obj.px_range = 2*pi*simgrid.symrange(obj.N(2))/(dx*obj.N(2));
            obj.py_range = 2*pi*simgrid.symrange(obj.N(1)).'/(dx*obj.N(1));
            obj.pz_range = reshape(2*pi*simgrid.symrange(obj.N(3))/(dx*obj.N(3)), [1, 1, obj.N(3)]);
       end
       function retval = p2(obj)
          %returns a 2 or 3 dimensional array with the magnitude of p squared.
          f_p2 = @(px, py) px.^2 + py.^2;
          retval = bsxfun(f_p2, obj.px_range, obj.py_range);
          if (obj.dimension == 3)
             f_pz = @(pxy2, pz) pxy2 + pz.^2;
             retval = bsxfun(f_pz, retval, obj.pz_range);
          end
       end
    end
	methods(Static)
		function range = symrange(N)
            if N==1
                range = 0;
            else
                range = fftshift(-N/2:N/2-1);
            end
        end
    end
end


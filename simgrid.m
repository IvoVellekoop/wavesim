classdef simgrid
    %grid parameters for simulation
    % Ivo M. Vellekoop
    properties
        N  % number of grid points in y,x-dimension
        dx  % grid resolution (x and y)
        padding % amount of added grid points on the x_boundary
        x_range % grid point coordinates in x-dimension
        y_range % grid point coordinates in x-dimension
        px_range % grid point coordinates in Fourier transformed x-dimension
        py_range % grid point coordinates in Fourier transformed y-dimension
    end
    methods
       function obj = simgrid(min_size, dx)
            % Construct a wave simulation grid object with specified width,
            % min_size = [height,width], minimum required size (will be rounded up)
			% dx = step size of grid (in arbitrary units)
			%% setup coordinates
            obj.N = pow2(nextpow2(min_size));
            obj.padding = obj.N - min_size; %total amoung of padding. Usually placed at right and bottom sides only (non-centric)
            obj.dx = dx;
            obj.x_range = (0:(obj.N(2)-1))*dx; %obj.x_range = linspace(min_x,max_x,obj.N(2));
            obj.y_range = (0:(obj.N(1)-1)).'*dx; %obj.y_range = linspace(min_y,max_y,obj.N(1));
            obj.px_range = 2*pi*simgrid.symrange(obj.N(2))/(dx*obj.N(2));
            obj.py_range = 2*pi*simgrid.symrange(obj.N(1)).'/(dx*obj.N(1));
       end
    end
	methods(Static)
		function range = symrange(N)
			range = fftshift(-N/2:N/2-1);
		end;
    end
end


classdef simgrid
    %grid parameters for simulation
    % Ivo M. Vellekoop
    properties
        Nx  % number of grid points in x-dimension
        Ny  % number of grid points in y-dimension
        dx  % grid resolution (x and y)
        x_margin % amount of grid points on the x_boundary
        y_margin % amount of grid points on the y_boundary
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
            obj.N = 2^nextpow2(min_size);
            obj.padding = obj.N - min_size; %total amoung of padding. Usually placed at right and bottom sides only (non-centric)
            obj.dx = dx;
            obj.x_range = simgrid.symrange(obj.N(2))*dx;
            obj.y_range = simgrid.symrange(obj.N(1)).'*dx;
            obj.px_range = 2*pi*simgrid.symrange(obj.N(2))/(dx*obj.Nx);
            obj.py_range = 2*pi*simgrid.symrange(obj.N(1)).'/(dx*obj.Ny);
       end
	methods(Static)
		function range = symrange(N)
			return fftshift(-N/2:N/2-1);
		end;
    end
end


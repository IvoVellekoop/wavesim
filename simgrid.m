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
       function obj = simgrid(Nx, Ny, dx)
            % Construct a wave simulation grid object with specified width,
            % height, and resolution. The wavelength is fixed at 1.
            %% setup coordinates
            obj.Nx = 2 ^ ceil(log2(Nx));
            obj.Ny = 2 ^ ceil(log2(Ny));
            
            obj.x_margin = (obj.Nx - Nx)/2;
            obj.y_margin = (obj.Ny - Ny)/2;
  
            obj.dx = dx;
            obj.x_range = fftshift(-obj.Nx/2:obj.Nx/2-1)*dx;
            obj.y_range = fftshift(-obj.Ny/2:obj.Ny/2-1).'*dx;
            obj.px_range = 2*pi*fftshift(-obj.Nx/2:obj.Nx/2-1)/(dx*Nx);
            obj.py_range = 2*pi*fftshift(-obj.Ny/2:obj.Ny/2-1).'/(dx*Ny);
            
            disp(['data size set to [',num2str(obj.Nx),',',num2str(obj.Ny),']']);
       end
    end
end


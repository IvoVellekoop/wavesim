classdef simgrid
    %grid parameters for simulation
    % Ivo M. Vellekoop
    properties
        Nx  % number of grid points in x-dimension
        Ny  % number of grid points in y-dimension
        dx  % grid resolution (x and y)
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
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.dx = dx;
            obj.x_range = fftshift(-Nx/2:Nx/2-1)*dx;
            obj.y_range = fftshift(-Ny/2:Ny/2-1).'*dx;
            obj.px_range = 2*pi*fftshift(-Nx/2:Nx/2-1)/(dx*Nx);
            obj.py_range = 2*pi*fftshift(-Ny/2:Ny/2-1).'/(dx*Ny);
       end
    end
end


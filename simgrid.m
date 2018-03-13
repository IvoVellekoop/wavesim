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
       function obj = simgrid(min_size, dx, pad)
            % Construct a wave simulation grid object with specified size (2-D or 3-D)
            % as a convention, the first, second, and third indices address the
            % x, y, and z dimensions, respectively.
            %
            % min_size = [height,width], minimum required size (will be rounded up)
			% dx = step size of grid (in arbitrary units)
            % pad = logic vector indicating whether or not to pad the
            % simulation size up to a size that is convenient for fft
			%% setup coordinates
            assert(isrow(min_size) && numel(min_size) == 3 && isrow(pad) && numel(pad) == 3);
            
            N = min_size;
            N(pad) = simgrid.efficient_size(min_size(pad));
            
            obj.padding = N - min_size; %total amount of zero padding. Placed at right and bottom sides only (non-centric)
            obj.dx = dx;
            obj.x_range = reshape((0:(N(1)-1))*dx, [N(1), 1, 1]);
            obj.y_range = reshape((0:(N(2)-1))*dx, [1, N(2), 1]);
            obj.z_range = reshape((0:(N(3)-1))*dx, [1, 1, N(3)]);
            obj.px_range = reshape(2*pi*simgrid.symrange(N(1))/(dx*N(1)), [N(1), 1, 1]);
            obj.py_range = reshape(2*pi*simgrid.symrange(N(2))/(dx*N(2)), [1, N(2), 1]);
            obj.pz_range = reshape(2*pi*simgrid.symrange(N(3))/(dx*N(3)), [1, 1, N(3)]);
            obj.N = N;
       end
    end
	methods(Static)
        function range = symrange(N)
            % constructs a range of numbers around 0 that is compatible
            % with the coordinates after a fft.
            % When N is odd, the range is symmetric around 0 (from -floor(N/2) to floor(N/2)).
            % When N is even, the range is almost symmetric and includes 0
            % (from -N/2 to N/2-1)
            % 
            range = fftshift((0:N-1)-floor(N/2));
        end
        function sz = efficient_size(min_size)
            % returns nearest size greater than or equal to min_size
            % for which the fft is efficient.
            % cuFFT is efficient for sizes that can be factored as 2^a * 3^b * 5^c * 7^d
            %
            sz = min_size;
            for s_i=1:length(sz)
                s = sz(s_i);
                f = factor(s);
                while f(end) > 7
                    s = s + 1;
                    f = factor(s);
                end
                sz(s_i) = s;
            end
        end
    end
end


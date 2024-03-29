classdef SimulationGrid
    % grid parameters for simulation
    % Ivo Vellekoop & Gerwin Osnabrugge 2016-2020
    properties
        N  % number of grid points in y,x,z-dimension
        dx  % grid resolution (same in all directions)
        padding % amount of added zero padding
        periodic % indicates whether simulation is periodic in y,x,z
        x_range % grid point coordinates in x-dimension
        y_range % grid point coordinates in y-dimension
        z_range % grid point coordinates in z-dimension
        dpx % step size in px
        dpy % step size in py
        dpz % step size in pz
        px_range % grid point coordinates in Fourier transformed x-dimension
        py_range % grid point coordinates in Fourier transformed y-dimension
        pz_range % grid point coordinates in Fourier transformed z-dimension
    end
    methods
       function obj = SimulationGrid(min_size, dx, periodic)
            % Construct a wave simulation grid object with specified size (2-D or 3-D)
            % as a convention, the first, second, and third indices address the
            % x, y, and z dimensions, respectively.
            %
            % min_size = [height,width], minimum required size (will be rounded up)
			% dx = step size of grid (in arbitrary units)
            % periodic = logic vector indicating whether or not to pad the
            % simulation size up to a size that is convenient for fft
			%% setup coordinates
            validateattributes(min_size, {'numeric'}, {'positive', 'size', [1,3]}); 
            validateattributes(dx, {'numeric'}, {'positive', 'scalar'});
            validateattributes(periodic, {'logical'}, {'size', [1,3]}); 
            
            min_size = ceil(min_size); %round up minimum required size to positive integer
            
            N = min_size;
            N(~periodic) = SimulationGrid.efficient_size(min_size(~periodic)); %increase size to efficient number for fft
            
            obj.periodic = periodic;
            obj.padding = N - min_size; %total amount of zero padding. Placed at right and bottom sides only (non-centric)
            obj.dx = dx;
            obj.x_range = reshape((0:(N(2)-1))*dx, [1, N(2), 1]);
            obj.y_range = reshape((0:(N(1)-1))*dx, [N(1), 1, 1]);
            obj.z_range = reshape((0:(N(3)-1))*dx, [1, 1, N(3)]);
            obj.dpx = 2*pi/(dx*N(2));
            obj.dpy = 2*pi/(dx*N(1));
            obj.dpz = 2*pi/(dx*N(3));
            obj.px_range = reshape(obj.dpx*SimulationGrid.symrange(N(2)), [1, N(2), 1]);
            obj.py_range = reshape(obj.dpy*SimulationGrid.symrange(N(1)), [N(1), 1, 1]);
            obj.pz_range = reshape(obj.dpz*SimulationGrid.symrange(N(3)), [1, 1, N(3)]);
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
            % cuFFT is efficient for sizes that can be factored as 2^a *
            % 3^b * 5^c * 7^d * 11^e (tested empirically using test_efficient_size)
            %
            sz = min_size;
            for s_i=1:length(sz)
                s = sz(s_i);
                f = factor(s);
                while f(end) > 11 || (length(f)>2 && f(end-1) > 5)
                    s = s + 1;
                    f = factor(s);
                end
                sz(s_i) = s;
            end
        end
        function test_efficient_size()
            s_range = 1:200;
            times = zeros(numel(s_range), 1);
            efficient = zeros(numel(s_range), 1, 'logical');
            for s_i=1:numel(s_range)
                s = s_range(s_i);
                f = factor(s);
                efficient(s_i) = true;
                if f(end) > 11 || (length(f)>2 && f(end-1) > 5)
                    efficient(s_i) = false;
                end
                data = gpuArray((1+1i)*ones([s,s,s], 'single'));
                times(s_i) = gputimeit(@() fftn(data)) / (s * log2(s));
                plot(s_range(efficient), times(efficient), '.'); 
                hold on;
                plot(s_range(~efficient), times(~efficient), '*'); 
                hold off;
                drawnow;
            end
        end
    end
end


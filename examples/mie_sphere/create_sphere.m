function [n, x_range, y_range, z_range] = create_sphere(options, N, oversample)
%% Returns a 3-D matrix of refractive indices for single sphere.
% The refractive index will be options.sphere_index inside a sphere
% and options.bg_index outside the sphere.
% Other parameters:
% options.sphere_radius The radius of the sphere (in micrometers)
% options.pixel_size    The size of a single pixel (in micrometers)
% N                     size of the matrix that is returned (in pixels) may be scalar
%                       or 3-element vector)
%
% returns coordinate vectors for x,y,z (in micrometers)

if numel(N)==1
    N = ones(1,3)*N;
end
if ~isfield(options, 'center')
    options.center = [0,0,0];
end
if nargin > 2
    opt = options;
    os_range = linspace(-0.5+0.5/oversample, 0.5-0.5/oversample, oversample);
    eps = zeros(N);
    for ox = os_range
        for oy = os_range
            for oz = os_range
                opt.center = options.center + [ox, oy, oz];
                eps = eps + create_sphere(opt, N).^2;
            end
        end
    end
    [~, x_range, y_range, z_range] = create_sphere(options, N);
    n = sqrt(eps / oversample^3);
    return;
end

%% Calculate coordinate ranges
center  = (N/2)*options.pixel_size + options.center;
x_range = reshape((1:N(2))*options.pixel_size - center(1), [1,N(2),1]);
y_range = reshape((1:N(1))*options.pixel_size - center(2), [N(1),1,1]);
z_range = reshape((1:N(3))*options.pixel_size - center(3), [1,1,N(3)]);

%% Calculate refractive index
inside = (x_range.^2 + y_range.^2 + z_range.^2) < options.sphere_radius.^2;
n = ones(N) * options.bg_index...
    + inside * (options.sphere_index - options.bg_index);



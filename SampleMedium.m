function obj = SampleMedium(refractive_index, options)
% SampleMedium - Generates a sample object for use in wave
% simulations (wavesim, PSTD, FDTD)
% internally, the object stores a map of the relative dielectric
% constant (e_r), which is the refractive index squared. The e_r map
% is padded with absorbing boundaries, and then expanded to the next
% multiple of 2^N in each dimension (for fast fourier transform)
%
% refractive_index   = refractive index map, may be complex, need not
%                      be square
% options.pixel_size = size of a grid pixel
% options.boundary_widths = vector with widths of the absorbing
%                          boundary (in pixels) for each dimension. Set element
%                          to 0 for periodic boundary.
% options.boundary_strength = maximum value of imaginary part of e_r
%                             in the boundary
%
% returned object:
% obj.e_r = relative dielectric constand map, with absorbing
%           boundaries appended, and padded to nearest multiple of 2^k,2^m
% obj.e_r_max = maximum real part of obj.e_r
% obj.e_r_min = minimum real part of obj.e_r
% obj.e_r_center = (e_r_max+e_r_min)/2
% obj.roi = size of the original refractive index map without padding
% obj.grid = simgrid object, with x and k ranges
% Saroch Leedumrongwatthanakun 2015

%calculate size of dielectric constant map
B = options.boundary_widths;
S = size(refractive_index);
obj.grid = simgrid(S+2*B, options.pixel_size); %padds to next power of 2 in each dimension
obj.roi = cell(2);
obj.roi{1} = B(1)+(1:S(1)); %the region of interest is between the boundaries
obj.roi{2} = B(2)+(1:S(2)); %the region of interest is between the boundaries

%calculate e_r and min/max values
obj.e_r = refractive_index.^2;
obj.e_r_min = min(real(obj.e_r(:)));
obj.e_r_max = max(real(obj.e_r(:)));
obj.e_r_center = (obj.e_r_min + obj.e_r_max)/2;

%add absorbing boundaries
obj.e_r = padarray(obj.e_r, B, 'replicate', 'both'); %absorbing boundaries
x = [(B(2):-1:1), zeros(1,S(2)), (1:B(2))] / max(B(2),1);
y = [(B(1):-1:1), zeros(1,S(1)), (1:B(1))].' / max(B(1),1);
BG = 1.0i*options.boundary_strength;
f_boundary = @(x, y) min(x.^2 + y.^2, 1).^2;
obj.e_r = obj.e_r + BG * bsxfun(f_boundary, x, y);

%add zero padding
obj.e_r = padarray(obj.e_r, obj.grid.N - size(obj.e_r), BG, 'post');
end

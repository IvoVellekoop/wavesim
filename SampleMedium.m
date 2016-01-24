function obj = SampleMedium(refractive_index, options)
% SampleMedium - Generates a sample object for use in wave
% simulations (wavesim, PSTD, FDTD)
% 
% internally, the object stores a map of the relative dielectric
% constant (e_r), which is the refractive index squared. The e_r map
% is padded with absorbing boundaries, and then expanded to the next
% multiple of 2^N in each dimension (for fast fourier transform)
%
% refractive_index   = refractive index map, may be complex, need not
%                      be square
% options.pixel_size = size of a grid pixel in any desired unit (e.g.
% microns)
% options.lambda     = wavelength (in same units as pixel_size)
% options.boundary_widths = vector with widths of the absorbing
%                          boundary (in pixels) for each dimension. Set element
%                          to 0 for periodic boundary.
% options.boundary_strength = maximum value of imaginary part of e_r
%                             in the boundary
% options.boundary_type = boundary type. Currently supports 'parabola' and
% 'PML' (default)
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

%% Set default values
if (~isfield(options, 'boundary_type'))
    options.boundary_type = 'PML2';
end;



%%calculate size of dielectric constant map
B = options.boundary_widths;
S = size(refractive_index);
obj.grid = simgrid(S+2*B, options.pixel_size); %padds to next power of 2 in each dimension
obj.roi = cell(3, 1);
obj.roi{1} = B(1)+(1:S(1)); %the region of interest is between the boundaries
obj.roi{2} = B(2)+(1:S(2)); %the region of interest is between the boundaries
if (length(S) == 3)
    obj.roi{3} = B(3)+(1:S(3)); %the region of interest is between the boundaries
else
    obj.roi{3} = 1;
end

%% Currently, the simulation will always be padded to a power of 2
% This is ok if we have boundaries, but if we have periodic boundary
% conditions, the field must have a size that is a power of 2 already
if (B(1)==0 && S(1) ~= pow2(nextpow2(S(1)))) || (B(2)==0 && S(2) ~= pow2(nextpow2(S(2))))
    error('If periodic boundary conditions are used, the sample size should be a power of 2 in that direction');
end

%calculate e_r and min/max values
obj.e_r = refractive_index.^2;
obj.e_r_min = min(real(obj.e_r(:)));
obj.e_r_max = max(real(obj.e_r(:)));
obj.e_r_center = (obj.e_r_min + obj.e_r_max)/2;

%add absorbing boundaries
%the shape of the boundary is determined by f_boundary_curve, a function
%that takes a position (in pixels, 0=start of boundary) and returns
%Delta e_r for the boundary. 
%
Bmax = max(B); %used to calculate expected amount of leakage through boundary
k0=0; c=0;
if (strcmp(options.boundary_type(1:3), 'PML')) %common stuff for all PMLs
    %perfectly matched layer boundary conditions.
    % the base refractive index of the boundaries is chosen as
    %the average value of e_r at all 4 boundaries
    %ignore boundaries with periodic boundary conditions
    e_sum = 0;
    e_count = 0;
    if (B(1)) 
        e_sum = e_sum + sum(obj.e_r(1,:)) + sum(obj.e_r(end,:));
        e_count = e_count + 2*size(obj.e_r, 2);
    else
        e_sum = e_sum + sum(obj.e_r(:,1)) + sum(obj.e_r(:,end));
        e_count = e_count + 2*size(obj.e_r, 1);
    end;
    e_0 = e_sum/e_count;
    obj.e_r = padarray(obj.e_r, B, e_0, 'both'); 
    k0 = sqrt(e_0)*2*pi/ (options.lambda / options.pixel_size); %k0 in 1/pixels
    % maximum value of the boundary (see Mathematica file = c(c-2ik0) = boundary_strength)
    % ||^2 = |c|^2 (|c|^2 + 4k0^2)   [assuming c=real, better possible when c complex?]
    % when |c| << |2k0| we can approximage: boundary_strength = 2 k0 c
    c = options.boundary_strength*k0^2 / (2*k0);
end
switch (options.boundary_type)
    %case 'PML' %Nth order smooth?
    %    N=3;
    %    f_boundary_curve = @(r) 1/k0^2*(c^(N+2)*r.^N.*(N+1.0+(2.0i*k0-c)*r)) ./ (factorial(N)*exp(c*r));
    %    obj.leakage = exp(-c*Bmax)*exp(c*Bmax);
    case 'PML5' %5th order smooth
        f_boundary_curve = @(r) 1/k0^2*(c^7*r.^5.*(6.0+(2.0i*k0-c)*r)) ./ (720+720*c*r+360*c^2*r.^2+120*c^3*r.^3+30*c^4*r.^4+6*c^5*r.^5+c^6*r.^6);
        obj.leakage = exp(-c*Bmax)*(720+720*c*Bmax+360*c^2*Bmax.^2+120*c^3*Bmax.^3+30*c^4*Bmax.^4+6*c^5*Bmax.^5+c^6*Bmax.^6)/24;
    case 'PML4' %4th order smooth
        f_boundary_curve = @(r) 1/k0^2*(c^6*r.^4.*(5.0+(2.0i*k0-c)*r)) ./ (120+120*c*r+60*c^2*r.^2+20*c^3*r.^3+5*c^4*r.^4+c^5*r.^5);
        obj.leakage = exp(-c*Bmax)*(120+120*c*Bmax+60*c^2*Bmax.^2+20*c^3*Bmax.^3+5*c^4*Bmax.^4+c^5*Bmax.^5)/24;
    case 'PML3' %3rd order smooth
        f_boundary_curve = @(r) 1/k0^2*(c^5*r.^3.*(4.0+(2.0i*k0-c)*r)) ./ (24+24*c*r+12*c^2*r.^2+4*c^3*r.^3+c^4*r.^4);
        obj.leakage = exp(-c*Bmax)*(24+24*c*Bmax+12*c^2*Bmax.^2+4*c^3*Bmax.^3+c^4*Bmax.^4)/24;
    case 'PML2' %2nd order smooth
        f_boundary_curve = @(r) 1/k0^2*(c^4*r.^2.*(3.0+(2.0i*k0-c)*r)) ./ (6+6*c*r+3*c^2*r.^2+c^3*r.^3);
        obj.leakage = exp(-c*Bmax)*(6+6*c*Bmax+3*c^2*Bmax.^2+c^3*Bmax.^3)/6;
    case 'PML1' %1st order smooth
        f_boundary_curve = @(r) 1/k0^2*(c^3*r.*(2.0+(2.0i*k0-c)*r)) ./ (2.0+2.0*c*r+c^2*r.^2) / k0^2; %(divide by k0^2 to get relative e_r)
        obj.leakage = exp(-c*Bmax)*(2+2*c*Bmax+c^2*Bmax.^2)/2;
    case 'parabola'
        %original boundary conditions.
        %the real part of e_r is replicated from the edge pixels of the
        %structure to reduce reflections at the edges
        %then, an imaginary part (absorption) is added.
        %The amount of absorption increases quadratically with depth into
        %the layer, reaching a maximum value of 'options.boundary_strength'
        %at the very edge
        obj.e_r = padarray(obj.e_r, B, 'replicate', 'both');
        R = max(max(B(2), B(1)),1);
        f_boundary_curve = @(r) 1.0i*options.boundary_strength * min(r/R, 1).^2;
    otherwise
        error(['unknown boundary type' obj.boundary_type]);
end;
    
%calculate boundaries (todo: implement for 3d)
if (sum(B) > 0)
    x = [(B(2):-1:1), zeros(1,S(2)), (1:B(2))];
    y = [(B(1):-1:1), zeros(1,S(1)), (1:B(1))].';
    if (length(B) == 3)
        z = [(B(3):-1:1), zeros(1,S(3)), (1:B(3))];
        z = reshape(z, [1,1,length(z)]);
        obj.e_r = obj.e_r + f_boundary_curve(sqrt(simgrid.dist2_3d(x,y,z)));
    else
        obj.e_r = obj.e_r + f_boundary_curve(sqrt(simgrid.dist2_2d(x,y)));
    end
end

%add padding. Note that we are padding with e_r_center. This is _not_
%optimal (optimal is to choose the boundaries to fill the full simulation
%grid). Padding with 0 may be better, but in that case the amount of
%boundary leakage strongly depends on the amount of padding (and hence, on
%the size of the simulation). This behavior is not desired at the moment.
%
obj.e_r = padarray(obj.e_r, obj.grid.Nred - size(obj.e_r), obj.e_r_center, 'post');
end

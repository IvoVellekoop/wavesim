%%% Script used to calculate the mie sphere analytical solution using Mie
%%% theory. May take several minutes to finish
%%% note: This script makes use of the external matlab code package 'MatScat',
%%% which can be found at: mathworks.com/matlabcentral/fileexchange/36831-matscat

%% Add pathway to MatScat code folder
dir_script = fileparts(mfilename('fullpath'));
if exist([dir_script,'\MatScat'],'dir') % check if MatScat package folder is correctly copied
    addpath(genpath('MatScat'));
else
    error(['MatScat package not included for copyright reasons.', newline, ...
           'This package can be found at: mathworks.com/matlabcentral/fileexchange/36831-matscat', newline, ...
           'Make sure to copy the MatScat package to the "wavesim/examples/mie_sphere/" folder']);
end

%% Sphere parameters               
opt.PPW = 5;                                        % number of grid points per wavelength
opt.lambda = 1;                                     % wavelength in vacuum (in um)
opt.sphere_radius = 6;                              % radius of the scattering sphere (in um)
opt.sphere_index = 1.2;                             % refractive index of the scattering sphere
opt.bg_index = 1.00;                                % refractive index of surrounding medium
opt.pixel_size = opt.lambda/opt.PPW;                % pixel spacing (in um)
N = [120,120,120];                                  % size of grid (in pixels)

%% Create Mie sphere refractive index distribution and corresponding grid
[n_sphere, x_r, y_r, z_r] = create_sphere(opt, N);
[xc, yc, zc] = ndgrid(x_r, y_r, z_r);

%% calculate analytical results (see MatScat folder for more details)
E_theory = calcmie_nf(opt.sphere_radius, opt.sphere_index, opt.bg_index, opt.lambda, xc, yc, zc);

% for some reason is the center pixel value for the returned field is zero
% when even number of grid points are used. we fill in center pixel with 
% the average value of direct neighbours
if mod(N(1),2) == 0
    E_theory(end/2,end/2,end/2,1) = (E_theory(end/2,end/2-1,end/2,1)+E_theory(end/2,end/2+1,end/2,1))/2;
    E_theory(end/2,end/2,end/2,2) = (E_theory(end/2,end/2-1,end/2,2)+E_theory(end/2,end/2+1,end/2,2))/2;
    E_theory(end/2,end/2,end/2,3) = (E_theory(end/2,end/2-1,end/2,3)+E_theory(end/2,end/2+1,end/2,3))/2;
end

% change order of dimensions to Matlab standard [y,x,z,pol]
E_theory = permute(E_theory,[2 1 3 4]);

%% save result
save('mie_sphere_theory.mat','E_theory','n_sphere','opt', 'x_r','y_r','z_r','N');

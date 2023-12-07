% Run tests for the calccyl_nf function
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Refresh workspace
close all;
clear variables;

%% Define test parameters
%dia = 2e-6;         % cylinder diameter
%ns = 1.33;          % cylinder refractive index (complex)
nm = 1.52;          % outer medium refractive index (real)
lambda = 600e-9;    % vacuum wavelength
nang = 1800;        % number of far field angles to evaluate
%zeta = 90;          % cylinder inclination angle (90 deg = perpendicular)

% stratified cylinder
dia = [1e-6, 2e-6]; % cylinder diameter
ns = [1.5, 1.33];   % cylinder refractive index (complex)
zeta = 90;          % cylinder inclination angle (90 deg = perpendicular)

sx = 2*dia(end);    % size of grid in x
sy = 2*dia(end);    % size of grid in y

Nx = 100;           % number of grid points in x
Ny = 100;           % number of grid points in y

conv = 1;           % convergence factor
tf_flag = false;    % total field flag
cc_flag = true;     % cartesian coordinates flag

%% Set up near field coordinates
deltax = sx/Nx;
deltay = sy/Ny;
nx = ((0:(Nx - 1)) - Nx/2.)*deltax;
ny = ((0:(Ny - 1)) - Ny/2.)*deltay;
[xf, yf] = ndgrid(nx, ny);
zf = zeros(size(xf));

%% Calculate near field solution
tic
[Ep, En, Hp, Hn, Sp, Sn, T, C, ang] = calccyl_nf(dia/2., ns, nm, ...
    lambda, xf, yf, zf, zeta, ...
    'ConvergenceFactor', conv, ...
    'TotalField', tf_flag, ...
    'Cartesian', cc_flag, ...
    'nang', nang);
toc

%% Plot near field solution
fields = {En(:,:,1), En(:,:,2), Ep(:,:,3), Hp(:,:,1), Hp(:,:,2), ...
    Hn(:,:,3)};

if cc_flag
    fldttl = {'E_x', 'E_y', 'E_z', 'H_x', 'H_y', 'H_z'};
else %if cc_flag
    fldttl = {'E_{rho}', 'E_{phi}', 'E_z', 'H_{rho}', 'H_{phi}', ...
        'H_z'}; %#ok<*UNRCH>
end %if cc_flag

figure
for ifld=1:length(fldttl)
    subplot(2, 3, ifld);
    if ~isempty(fields{ifld})
        imagesc(flipud(rot90(abs(fields{ifld}).^2)));
    end %if ~isempty(fields{ifld})
    title(fldttl{ifld});
end %for ifld=1:length(fldlst)

%% Plot Poynting vector
figure()
subplot(1,2,1)
quiver(xf,yf,Sp(:,:,1), Sp(:,:,2))
axis image
subplot(1,2,2)
quiver(xf,yf,Sn(:,:,1), Sn(:,:,2))
axis image

%% Plot far field solution

% Differential scattering cross sections
fctr = 2/pi/C.k;
dCsdOp = fctr*squeeze(abs(T(1,1,:).^2));    % TM
dCsdOn = fctr*squeeze(abs(T(2,2,:).^2));    % TE
dCsdO = 0.5*(dCsdOp + dCsdOn);              % unpolarized

figure
subplot(1,3,1);
semilogy(ang, dCsdOp);
title('TM')
subplot(1,3,2);
semilogy(ang, dCsdOn);
title('TE')
subplot(1,3,3);
semilogy(ang, dCsdO);
title('unpolarised')

for i=1:3
    subplot(1,3,i);
    xlabel('Scattering angle [^\circ]')
    ylabel('Differential scattering cross section [m]')
    xlim([ang(1), ang(end)])
end %for i=1:3

disp('Cross sections:');
disp(C);
disp('Efficiencies:');
disp(getEfficiencies(C, dia(end)/2., 2));
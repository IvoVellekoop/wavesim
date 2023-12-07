% Run tests for the calcmie_nf function
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Refresh workspace
close all;
clear variables;

%% Define test parameters
% dia = 2e-7;         % sphere diameter
% ns = 1.33;          % sphere refractive index (complex)

% stratified sphere
dia = [1e-7, 2e-7, 3e-7, 4e-7];   % sphere diameter
ns = [1.1, 1.33, 1.8, 1.33];      % sphere refractive index (complex)

dia = [2e-7, 3e-7, 4e-7]; % sphere diameter
ns = [1., 1.3, 1.8];      % sphere refractive index (complex)

nm = 1.;            % outer medium refractive index (real)
lambda = 600e-9;    % vacuum wavelength
nang = 1800;        % number of far field angles to evaluate

sx = 2*dia(end);    % size of grid in x
sy = 2*dia(end);    % size of grid in y

Nx = 100;           % number of grid points in x
Ny = 100;           % number of grid points in y

conv = 1;           % convergence factor
tf_flag = true;     % total field flag
cc_flag = true;     % cartesian coordinates flag

rad = dia/2.;       % sphere radius

%% Set up near field coordinates
deltax = sx/Nx;
deltay = sy/Ny;
nx = ((0:(Nx - 1)) - Nx/2.)*deltax;
ny = ((0:(Ny - 1)) - Ny/2.)*deltay;
[xf, yf] = ndgrid(nx, ny);
zf = zeros(size(xf));

%% Calculate near field solution
tic
[E, H, P, S, C, ang] = calcmie_nf(dia/2., ns, nm, ...
    lambda, xf, yf, zf, ...
    'ConvergenceFactor', conv, ...
    'TotalField', tf_flag, ...
    'Cartesian', cc_flag, ...
    'nang', nang);
toc

%% Plot near field solution
fields = {E(:,:,1), E(:,:,2), E(:,:,3), H(:,:,1), H(:,:,2), ...
    H(:,:,3)};

if cc_flag
    fldttl = {'E_x', 'E_y', 'E_z', 'H_x', 'H_y', 'H_z'};
else %if cc_flag
    fldttl = {'E_{rho}', 'E_{phi}', 'E_{theta}', 'H_{rho}', 'H_{phi}', ...
        'H_{theta}'}; %#ok<*UNRCH>
end %if cc_flag

figure
for ifld=1:length(fldttl)
    subplot(2, 3, ifld);
    if ~isempty(fields{ifld})
        imagesc(nx, ny, flipud(rot90(abs(fields{ifld}).^2)));
        rectangle(...
            'Position', [-rad(end), -rad(end), dia(end), dia(end)], ...
            'Curvature', [1,1])
    end %if ~isempty(fields{ifld})
    title(fldttl{ifld});
end %for ifld=1:length(fldlst)

%% Plot Poynting vector
figure()
imagesc(nx, ny, sqrt(P(:,:,1).^2 + P(:,:,2).^2 + P(:,:,3).^2))
rectangle('Position', [-rad(end), -rad(end), dia(end), dia(end)], ...
    'Curvature', [1,1])
axis image

%% Plot far field solution

% Differential scattering cross sections
fctr = 2/pi/C.k;
dCsdOp = fctr*squeeze(abs(S(1,1,:).^2));    % parallel
dCsdOn = fctr*squeeze(abs(S(2,2,:).^2));    % perpendicular
dCsdO = 0.5*(dCsdOp + dCsdOn);              % unpolarized

figure
subplot(1,3,1);
semilogy(ang, dCsdOp);
title('parallel')
subplot(1,3,2);
semilogy(ang, dCsdOn);
title('perpendicular')
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
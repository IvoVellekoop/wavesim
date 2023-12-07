% Run tests for the calccyl_nf_multi function
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Refresh workspace
close all;
clear variables;

%% initialization values of the scatterer
dia = 2e-7;         % the cylinder diameter
ns = 1.33 +0.1j;    % the refractive index of the cylinder
nm = 1.;            % the refractive index of the outer medium
lambda = 600e-9;    % the wavelength
nang = 1800;        % number of far field angles to evaluate
iphi = 45;          % the incident angle

conv = 1;           % convergence factor
tf_flag = false;    % total field flag
cc_flag = true;     % cartesian coordinates flag

if cc_flag
    fldttl = {'E_x', 'E_y', 'E_z', 'H_x', 'H_y', 'H_z'};
else %if cc_flag
    fldttl = {'E_{rho}', 'E_{phi}', 'E_z', 'H_{rho}', 'H_{phi}', 'H_z'};
    %#ok<*UNRCH>
end %if cc_flag



%% Calculate single scattering

% Set up near field coordinates

orgn = [0., 0.];    % the center position of the scatterer

sx = 2*dia;         % size of grid in x
sy = 2*dia;         % size of grid in y

Nx = 100;           % number of grid points in x
Ny = 100;           % number of grid points in y

deltax = sx/Nx;
deltay = sy/Ny;
nx = ((0:(Nx - 1)) - Nx/2.)*deltax;
ny = ((0:(Ny - 1)) - Ny/2.)*deltay;
[xf, yf] = ndgrid(nx, ny);
zf = zeros(size(xf));

% Calculate single cylinder near field solution
[ Ep1, En1, Hp1, Hn1, Sp1, Sn1, T1, C1, ang1] = calccyl_nf( dia/2., ...
    ns, nm, lambda, xf-orgn(1), yf-orgn(2), zf, 90, ...
    'TotalField', tf_flag, ...
    'Cartesian', cc_flag, ...
    'nang', nang);

field1 = {En1(:,:,1), En1(:,:,2), Ep1(:,:,3), Hp1(:,:,1), Hp1(:,:,2), ...
    Hn1(:,:,3)};

% Plot single cylinder near field solution
figure
for ifld=1:length(fldttl)
    subplot(2, 3, ifld);
    if ~isempty(field1{ifld})
        imagesc(abs(flipud(rot90(field1{ifld}).^2)));
    end %if ~isempty(field1{ifld})
    title(fldttl{ifld});
end %for ifld=1:length(fldttl)

% Plot Poynting vector
figure
subplot(1,2,1)
quiver(xf, yf, Sp1(:,:,1), Sp1(:,:,2))
axis image
subplot(1,2,2)
quiver(xf, yf, Sn1(:,:,1), Sn1(:,:,2))
axis image

% Calculate multiple cylinder near field solution
[ Ep2, En2, Hp2, Hn2, Sp2, Sn2, T2, C2, ang2 ] = calccyl_nf_multi( ...
    dia/2., ns, nm, lambda, orgn(1), orgn(2), xf, yf, ...
    'iphi', iphi, ...
    'TotalField', tf_flag, ...
    'Cartesian', cc_flag, ...
    'nang', nang);
field2 = {En2(:,:,1), En2(:,:,2), Ep2(:,:,3), Hp2(:,:,1), Hp2(:,:,2), ...
    Hn2(:,:,3)};

% Plot multiple cylinder near field solution
figure
for ifld=1:length(fldttl)
    subplot(2, 3, ifld);
    if ~isempty(field2{ifld})
        imagesc(abs(flipud(rot90(field2{ifld}).^2)));
    end %if ~isempty(field2{ifld})
    title(fldttl{ifld});
end %for ifld=1:length(fldttl)

% Plot Poynting vector
figure
subplot(1,2,1)
quiver(xf, yf, Sp2(:,:,1), Sp2(:,:,2))
axis image
subplot(1,2,2)
quiver(xf, yf, Sn2(:,:,1), Sn2(:,:,2))
axis image

%% Plot far field solution

% Differential scattering cross sections
fctr1 = 2/pi/C1.k;
dCsdOp1 = fctr1*squeeze(abs(T1(1,1,:).^2));    % TM
dCsdOn1 = fctr1*squeeze(abs(T1(2,2,:).^2));    % TE
dCsdO1 = 0.5*(dCsdOp1 + dCsdOn1);              % unpolarized

figure
subplot(1,3,1);
semilogy(ang1, dCsdOp1);
title('TM')
subplot(1,3,2);
semilogy(ang1, dCsdOn1);
title('TE')
subplot(1,3,3);
semilogy(ang1, dCsdO1);
title('unpolarised')

for i=1:3
    subplot(1,3,i);
    xlabel('Scattering angle [^\circ]')
    ylabel('Differential scattering cross section [m]')
    xlim([ang2(1), ang2(end)])
end %for i=1:3

disp('Cross sections:');
disp(C1);
disp('Efficiencies:');
disp(getEfficiencies(C1, dia/2., 2));

% Differential scattering cross sections
fctr2 = 2/pi/C2.k;
dCsdOp2 = fctr2*squeeze(abs(T2(1,1,:).^2));    % TM
dCsdOn2 = fctr2*squeeze(abs(T2(2,2,:).^2));    % TE
dCsdO2 = 0.5*(dCsdOp2 + dCsdOn2);              % unpolarized

figure
subplot(1,3,1);
semilogy(ang2, dCsdOp2);
title('TM')
subplot(1,3,2);
semilogy(ang2, dCsdOn2);
title('TE')
subplot(1,3,3);
semilogy(ang2, dCsdO2);
title('unpolarised')

for i=1:3
    subplot(1,3,i);
    xlabel('Scattering angle [^\circ]')
    ylabel('Differential scattering cross section [m]')
    xlim([ang2(1), ang2(end)])
end %for i=1:3

disp('Cross sections:');
disp(C2);
disp('Efficiencies:');
disp(getEfficiencies(C2, dia/2., 2));

%% Calculate multiple cylinder scattering

% Set up near field coordinates

sx = 4*dia;         % size of grid in x
sy = 4*dia;         % size of grid in y

Nx = 100;           % number of grid points in x
Ny = 100;           % number of grid points in y

deltax = sx/Nx;
deltay = sy/Ny;
nx = ((0:(Nx - 1)) - Nx/2.)*deltax;
ny = ((0:(Ny - 1)) - Ny/2.)*deltay;
[xf, yf] = ndgrid(nx, ny);
zf = zeros(size(xf));

xc = [0., -1.5, 1.5]*dia;
yc = [0., .5, 0.]*dia;

% Calculate multiple cylinder near field solution
[ Ep1, En1, Hp1, Hn1, Sp1, Sn1, T1, C1, ang1 ] = calccyl_nf_multi( ...
    dia/2., ns, nm, lambda, xc, yc, xf, yf, ...
    'TotalField', tf_flag, ...
    'Cartesian', cc_flag, ...
    'nang', nang);

field1 = {En1(:,:,1), En1(:,:,2), Ep1(:,:,3), Hp1(:,:,1), Hp1(:,:,2), ...
    Hn1(:,:,3)};

% Plot single cylinder near field solution
figure
for ifld=1:length(fldttl)
    subplot(2, 3, ifld);
    if ~isempty(field1{ifld})
        imagesc(abs(flipud(rot90(field1{ifld}).^2)));
    end %if ~isempty(field1{ifld})
    title(fldttl{ifld});
end %for ifld=1:length(fldttl)

% Plot Poynting vector
figure
subplot(1,2,1)
quiver(xf, yf, Sp1(:,:,1), Sp1(:,:,2))
axis image
subplot(1,2,2)
quiver(xf, yf, Sn1(:,:,1), Sn1(:,:,2))
axis image


% Calculate multiple cylinder near field solution
[ Ep2, En2, Hp2, Hn2, Sp2, Sn2, T2, C2, ang2 ] = calccyl_nf_multi( ...
    dia/2., ns, nm, lambda, xc, yc, xf, yf, ...
    'iphi', iphi, ...
    'TotalField', tf_flag, ...
    'Cartesian', cc_flag, ...
    'nang', nang);
field2 = {En2(:,:,1), En2(:,:,2), Ep2(:,:,3), Hp2(:,:,1), Hp2(:,:,2), ...
    Hn2(:,:,3)};

% Plot multiple cylinder near field solution
figure
for ifld=1:length(fldttl)
    subplot(2, 3, ifld);
    if ~isempty(field2{ifld})
        imagesc(abs(flipud(rot90(field2{ifld}).^2)));
    end %if ~isempty(field2{ifld})
    title(fldttl{ifld});
end %for ifld=1:length(fldttl)

% Plot Poynting vector
figure
subplot(1,2,1)
quiver(xf, yf, Sp2(:,:,1), Sp2(:,:,2))
axis image
subplot(1,2,2)
quiver(xf, yf, Sn2(:,:,1), Sn2(:,:,2))
axis image

%% Plot far field solution

% Differential scattering cross sections
fctr1 = 2/pi/C1.k;
dCsdOp1 = fctr1*squeeze(abs(T1(1,1,:).^2));    % TM
dCsdOn1 = fctr1*squeeze(abs(T1(2,2,:).^2));    % TE
dCsdO1 = 0.5*(dCsdOp1 + dCsdOn1);              % unpolarized

figure
subplot(1,3,1);
semilogy(ang1, dCsdOp1);
title('TM')
subplot(1,3,2);
semilogy(ang1, dCsdOn1);
title('TE')
subplot(1,3,3);
semilogy(ang1, dCsdO1);
title('unpolarised')

for i=1:3
    subplot(1,3,i);
    xlabel('Scattering angle [^\circ]')
    ylabel('Differential scattering cross section [m]')
    xlim([ang2(1), ang2(end)])
end %for i=1:3

disp('Cross sections:');
disp(C1);
disp('Efficiencies:');
disp(getEfficiencies(C1, dia/2., 2));

% Differential scattering cross sections
fctr2 = 2/pi/C2.k;
dCsdOp2 = fctr2*squeeze(abs(T2(1,1,:).^2));    % TM
dCsdOn2 = fctr2*squeeze(abs(T2(2,2,:).^2));    % TE
dCsdO2 = 0.5*(dCsdOp2 + dCsdOn2);              % unpolarized

figure
subplot(1,3,1);
semilogy(ang2, dCsdOp2);
title('TM')
subplot(1,3,2);
semilogy(ang2, dCsdOn2);
title('TE')
subplot(1,3,3);
semilogy(ang2, dCsdO2);
title('unpolarised')

for i=1:3
    subplot(1,3,i);
    xlabel('Scattering angle [^\circ]')
    ylabel('Differential scattering cross section [m]')
    xlim([ang2(1), ang2(end)])
end %for i=1:3

disp('Cross sections:');
disp(C2);
disp('Efficiencies:');
disp(getEfficiencies(C2, dia/2., 2));


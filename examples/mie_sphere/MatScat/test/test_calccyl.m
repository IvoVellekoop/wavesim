% Run tests for the calccyl function
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Refresh workspace
close all;
clear variables;

%% Define test parameters
dia = 2*0.525e-6;  % cylinder diameter
ns = 1. + 0.j;     % cylinder refractive index (complex)
nm = 1.33;         % outer medium refractive index (real)
lambda = 632.8e-9; % vacuum wavelength
nang = 21;         % number of far field angles to evaluate
zeta = 90;         % cylinder inclination angle (90 deg = perpendicular)

% stratified cylinder
% dia = [1e-6, 2e-6]; % cylinder diameter
% ns = [1.474, 1.2];  % cylinder refractive index (complex)
% nm = 1.;            % outer medium refractive index (real)
% lambda = 600e-9;    % vacuum wavelength
% nang = 1800;        % number of far field angles to evaluate
% zeta = 90;          % cylinder inclination angle (90 deg = perpendicular)

%% Caclulate single cylinder far field solution
tic
[T, C, ang] = calccyl(dia/2., ns, nm, lambda, nang, zeta);
toc

%% Differential scattering cross sections
fctr = 2/pi/C.k;
dCsdOp = fctr*squeeze(abs(T(1,1,:).^2));    % TM
dCsdOn = fctr*squeeze(abs(T(2,2,:).^2));    % TE
dCsdO = 0.5*(dCsdOp + dCsdOn);              % unpolarized

%% Plot single cylinder far field solution
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
disp(getEfficiencies(C, dia/2., 2));

%% Calculate and plot Mueller matrix
M = getMuellerMatrix(T);
figure
idx = 1;
for i=1:4
    for j=1:4
        subplot(4,4,idx);
        if idx == 1
            semilogy(ang, squeeze(M(1,1,:)));
            ylabel('M_{11}')
        else %if idx == 1
            plot(ang, squeeze(M(i,j,:)./M(1,1,:)));
            ylabel(['M_{' num2str(i) num2str(j) '}/M_{11}'])
            ylim([-1, 1])
        end %if idx == 1
        xlabel('Scattering angle [^\circ]')
        xlim([ang(1), ang(end)])
        idx = idx + 1;
    end %for j=1:4
end %for i=1:4

T11 = squeeze(M(1,1,:)/M(1,1,1));
T33 = squeeze(M(3,3,:)./M(1,1,:));
T34 = squeeze(M(3,4,:)./M(1,1,:));
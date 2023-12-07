% Run tests for the calcmie function
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Refresh workspace
clear variables
close all;

%% Define test parameters
dia = 2e-6;         % sphere diameter
ns = 1.33 + .01j;   % sphere refractive index (complex)

% stratified sphere
% dia = [1e-6, 2e-6]; % sphere diameter
% ns = [1.33, 1.];    % sphere refractive index (complex)

nm = 1.52;          % outer medium refractive index (real)
lambda = 600e-9;    % vacuum wavelength

nang = 1800;        % number of far field angles to evaluate
conv = 1;           % convergence factor

rad = dia/2.;           % sphere radius
k = 2*pi/lambda*nm;    % wavenumber in medium n_m

%% Calculate amplitude scattering matrix
[S, C, ang] = calcmie(rad, ns, nm, lambda, nang, ...
    'ConvergenceFactor', conv);

%% Differential scattering cross sections
fctr = 2/pi/C.k;
dCsdOp = fctr*squeeze(abs(S(1,1,:).^2));    % parallel
dCsdOn = fctr*squeeze(abs(S(2,2,:).^2));    % perpendicular
dCsdO = 0.5*(dCsdOp + dCsdOn);              % unpolarized

%% Plot single cylinder far field solution
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


%% Calculate cross sections and efficiencies
Q = getEfficiencies(C, rad(end), 3);
disp('Cross sections:')
disp(C);
disp('Efficiencies:')
disp(Q);

%% Calculate and plot Mueller matrix
M = getMuellerMatrix(S);
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
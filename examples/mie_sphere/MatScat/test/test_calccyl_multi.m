% Run tests for the calccyl_multi function
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Refresh workspace
clear variables
close all;

%% Define test parameters
dia = 1e-6;         % cylinder diameter
ns = 1.33 + 1j;     % cylinder refractive index (complex)
nm = 1.52;          % outer medium refractive index (real)
lambda = 600e-9;    % vacuum wavelength
nang = 1800;        % number of far field angles to evaluate
iphi = 0;           % incident angle (degree)

%% Single scattering test
pos = [0, 0];       % cylinder positions

% Calculate with single cylinder function 
tic
[ Ts, Cs, angs ] = calccyl( dia/2., ns, nm, lambda, nang);
toc
fctrs = 2/pi/Cs.k;
dCsdOps = fctrs*squeeze(abs(Ts(1,1,:).^2));
dCsdOns = fctrs*squeeze(abs(Ts(2,2,:).^2));
dCsdOs = 0.5*(dCsdOps + dCsdOns);

% Calculate with multiple cylinder function
tic
[ Tm, Cm, angm ] = calccyl_multi( dia/2., ns, nm, lambda, ...
    pos(:,1), pos(:,2), nang, iphi);
toc
fctrm = 2/pi/Cm.k;
dCsdOpm = fctrm*squeeze(abs(Tm(1,1,:).^2));
dCsdOnm = fctrm*squeeze(abs(Tm(2,2,:).^2));
dCsdOm = 0.5*(dCsdOpm + dCsdOnm);

% Plot results
figure()
subplot(1,3,1);
semilogy(angs, dCsdOps);
hold on;
semilogy(angm, dCsdOpm, 'r');
title('TM')
subplot(1,3,2);
semilogy(angs, dCsdOns);
hold on;
semilogy(angm, dCsdOnm, 'r');
title('TE')
subplot(1,3,3);
semilogy(angs, dCsdOs);
hold on;
semilogy(angm, dCsdOm, 'r');
title('unpolarised')

for i=1:3
    subplot(1,3,i);
    xlabel('Scattering angle [^\circ]')
    ylabel('Differential scattering cross section [m]')
    xlim([angm(1), angm(end)])
end %for i=1:3

disp('Cross sections:');
disp(Cs);
disp(Cm);
disp('Efficiencies:');
disp(getEfficiencies(Cs, dia/2., 2));
disp(getEfficiencies(Cm, dia/2., 2));

%% Multiple Scattering test
pos = [0, 0; 1, 1; 2, 0]*dia;   % cylinder positions
dphi = 180;                     % incident angle shift

% Calculate with multiple cylinder function
tic
[ T1, C1, ang1 ] = calccyl_multi( dia/2., ns, nm, lambda, ...
    pos(:,1), pos(:,2), nang, iphi);
toc
fctr1 = 2/pi/C1.k;
dCsdOp1 = fctr1*squeeze(abs(T1(1,1,:).^2));
dCsdOn1 = fctr1*squeeze(abs(T1(2,2,:).^2));
dCsdO1 = 0.5*(dCsdOp1 + dCsdOn1);

% Calculate with multiple cylinder function and angle shift
tic
[ T2, C2, ang2 ] = calccyl_multi( dia/2., ns, nm, lambda, ...  
    pos(:,1), pos(:,2), nang, iphi+dphi);
toc
fctr2 = 2/pi/C2.k;
dCsdOp2 = fctr2*squeeze(abs(T2(1,1,:).^2));
dCsdOn2 = fctr2*squeeze(abs(T2(2,2,:).^2));
dCsdO2 = 0.5*(dCsdOp2 + dCsdOn2);

% Calculate with multiple cylinder function using Hankel function H_n^(2)
tic
[ T3, C3, ang3 ] = calccyl_multi( dia/2., conj(ns), nm, lambda, ...  
    pos(:,1), pos(:,2), nang, iphi, 'HankelFunction', 2);
toc
fctr3 = 2/pi/C3.k;
dCsdOp3 = fctr3*squeeze(abs(T3(1,1,:).^2));
dCsdOn3 = fctr3*squeeze(abs(T3(2,2,:).^2));
dCsdO3 = 0.5*(dCsdOp3 + dCsdOn3);

% Calculate with multiple cylinder function neglecting dependent scattering
tic
[ T4, C4, ang4 ] = calccyl_multi( dia/2., conj(ns), nm, lambda, ...  
    pos(:,1), pos(:,2), nang, iphi, 'DependentScattering', false);
toc

fctr4 = 2/pi/C4.k;
dCsdOp4 = fctr4*squeeze(abs(T4(1,1,:).^2));
dCsdOn4 = fctr4*squeeze(abs(T4(2,2,:).^2));
dCsdO4 = 0.5*(dCsdOp4 + dCsdOn4);

% Plot results

figure()
subplot(1,3,1);
semilogy(ang1, dCsdOp1);
hold on;
semilogy(ang2, dCsdOp2, 'r');
semilogy(ang3, dCsdOp3, 'g');
semilogy(ang4, dCsdOp4, 'm');
title('TM')
subplot(1,3,2);
semilogy(ang1, dCsdOn1);
hold on;
semilogy(ang2, dCsdOn2, 'r');
semilogy(ang3, dCsdOn3, 'g');
semilogy(ang4, dCsdOn4, 'm');
title('TE')
subplot(1,3,3);
semilogy(ang1, dCsdO1);
hold on;
semilogy(ang2, dCsdO2, 'r');
semilogy(ang3, dCsdO3, 'g');
semilogy(ang4, dCsdO4, 'm');
title('unpolarised')

for i=1:3
    subplot(1,3,i);
    xlabel('Scattering angle [^\circ]')
    ylabel('Differential scattering cross section [m]')
    xlim([angm(1), angm(end)])
end %for i=1:3
disp('Cross sections:');
disp(C1);
disp(C2);
disp(C3);
disp(C4);
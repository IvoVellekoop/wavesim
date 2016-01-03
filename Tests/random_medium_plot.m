addpath('..');
clear all; close all;
%% Load random medium simulation data
load('random_medium_results.mat');

%% plot results
fig_size = [8 7];
loglog(iterations_per_wavelength(2:end),errors_PSTD,'r+', 'MarkerSize', 10);
fixplot('Iterations per wavelength', 'Relative error', fig_size, '');

%% refractive index map image (real & imaginary)
x = 1:N(2)/PPW;
y = 1:N(1)/PPW;

figure(2);
imagesc(x,y,real(n_sample))
h = colorbar;
set(get(h,'Title'),'String','Re\{V\}','FontSize',16);
fixplot('x (\lambda)','y (\lambda)',fig_size,'');
axis square;

figure(3);
imagesc(x,y,imag(n_sample))
h = colorbar;
set(get(h,'Title'),'String','Im\{V\}','FontSize',16);
fixplot('x (\lambda)','y (\lambda)',fig_size,'');
axis square;

%% Field solutions (wavesim & PSTD)
figure(4);
imagesc(x,y,log(abs(E_wavesim)))
h = colorbar;
set(get(h,'Title'),'String','log|\psi|^2','FontSize',16);
fixplot('x (\lambda)','y (\lambda)',fig_size,'');
axis square;

figure(5);
imagesc(x,y,log(abs(E_PSTD{end})))
h = colorbar;
set(get(h,'Title'),'String','log|\psi|^2','FontSize',16);
fixplot('x (\lambda)','y (\lambda)',fig_size,'');
axis square;
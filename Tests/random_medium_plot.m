%%% script used for plotting the results of the random medium
%%% simulations.

close all;
addpath('..');
fp = '../../wavesimpaper/figures/'; % filepath for figures
%% Load random medium simulation data
load('random_medium_results.mat');

%% error
for t_i=2:length(dt_relative_range)
  errors_PSTD(t_i-1) = mean2(abs(E_PSTD{t_i-1} - E_wavesim).^2) / mean2(abs(E_wavesim).^2);
end

%% plot results
fig_size = [8 7];
loglog(iterations_per_wavelength(2:end),errors_PSTD,'r+', 'MarkerSize', 10);
fixplot('Iterations per wavelength', 'Relative error', fig_size, '');

print([fp 'random_results.eps'], '-depsc2'); % print figure to eps file
%% refractive index map image (real & imaginary)
x = 1:N(2)/PPW;
y = 1:N(1)/PPW;

% figure(2);
% imagesc(x,y,real(n_sample))
% h = colorbar;
% set(get(h,'Title'),'String','Re\{\epsilon_{r}\}','FontSize',16);
% fixplot('x (\lambda)','y (\lambda)',fig_size,'');
% axis square;
% colormap(jet);
% 
% figure(3);
% imagesc(x,y,imag(n_sample))
% h = colorbar;
% set(get(h,'Title'),'String','Im\{\epsilon_{r}\}','FontSize',16);
% fixplot('x (\lambda)','y (\lambda)',fig_size,'');
% axis square;
% colormap(jet);

%% Field solutions (wavesim & PSTD)
% wavesim solution
figure(4);
imagesc(x,y,log(abs(E_wavesim)))
h = colorbar;
set(get(h,'Title'),'String','log|\psi|^2','FontSize',16);
fixplot('x (\lambda)','y (\lambda)',fig_size,'');
axis square;
colormap(jet);

print([fp 'random_solution.eps'], '-depsc2'); % print figure to eps file

% PSTD solution
% figure(5);
% imagesc(x,y,log(abs(E_PSTD{end})))
% h = colorbar;
% set(get(h,'Title'),'String','log|\psi|^2','FontSize',16);
% fixplot('x (\lambda)','y (\lambda)',fig_size,'');
% axis square;
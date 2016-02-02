%%% script used for plotting the results of the random medium
%%% simulations.

close all;
addpath('..');
fp = '../../wavesimpaper/figures/'; % filepath for figures
%% Load random medium simulation data
load('disordered_medium_3d_results_subset.mat');

%% plot results
fig_size = [8 7];
theory = 1./iterations_per_wavelength.^4 * errors_PSTD(5) * iterations_per_wavelength(6).^4;
loglog(iterations_per_wavelength(2:end), theory(2:end), 'k', 'LineWidth', 2.0);
hold on
loglog(iterations_per_wavelength(2:end),errors_PSTD,'r+', 'MarkerSize', 10, 'LineWidth', 2.0);
fixplot('Iterations per wavelength', 'Relative error', fig_size, '');
set(gca,'YMinorGrid','on');
print([fp 'disordered_3d_results.eps'], '-depsc2'); % print figure to eps file

%% Field solutions (wavesim & PSTD)
% wavesim solution
fig_size = [8 7];
figure(4);
E = {E64, E48, E32, E1};
PPW = 4;
x = 1:128/PPW;
y = x;
%fixplot('', '', fig_size, '');
z = [64 48 32 1];
for p=1:4
    figure;
    imagesc(x,y,log(abs(E{p})))
    h = colorbar;
    set(get(h,'Title'),'String','log|\psi|^2','FontSize',16);
    fixplot('x (\lambda)','y (\lambda)',fig_size*0.75,'');
    axis square;
    colormap(jet);
    print([fp 'disordered_3d_solution_' num2str(z(p)) '.eps'], '-depsc2'); % print figure to eps file
end

% PSTD solution
% figure(5);
% imagesc(x,y,log(abs(E_PSTD{end})))
% h = colorbar;
% set(get(h,'Title'),'String','log|\psi|^2','FontSize',16);
% fixplot('x (\lambda)','y (\lambda)',fig_size,'');
% axis square;
%%% Script used to plot the accuracy of PSTD and wavesim in homogeneous
%%% medium
%Note: uncomment print commands to save figures in eps format
close all;
%% load figure_accuracy simulation results
load('figure_accuracy.mat');

%% Plot results
figure(1);

loglog(iterations_per_wavelength(1), relative_error(1), 's', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'k', 'MarkerSize', 8); % wavesim results
hold on;
loglog(iterations_per_wavelength(2:end), relative_error(2:end), 'r+', 'MarkerSize', 10, 'LineWidth', 2.0); % PSTD results
 
% plot line indicating expected error of PSTD
theory = 1./iterations_per_wavelength.^4 * relative_error(6) * iterations_per_wavelength(6).^4; 
loglog(iterations_per_wavelength(2:end), theory(2:end), 'k', 'LineWidth', 2.0);
loglog([0.1, 10E4], relative_error(1)*[1,1], '--','Color',[0,0.4,0.6]);
set(gca,'YMinorGrid','on');
hold off;

fixplot('Iterations per wavelength', 'Relative error', [8 7], '');
xlim([0.1, 1E5]);
ylim([1E-12, 10]);
h_legend = legend('Modified Born','PSTD');
set(h_legend,'FontSize',22);
fp = '../../figures/';
%print([fp 'accuracy_homogeneous.eps'], '-depsc2');
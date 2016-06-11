addpath('..');
load('homogeneous_medium_3d_results.mat');
thit = [iterations_per_wavelength(2)/2 iterations_per_wavelength(2:end-3)];
theory = 1./thit.^4 * relative_error(6) * iterations_per_wavelength(6).^4;

loglog(iterations_per_wavelength(1), relative_error(1), 's', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'k', 'MarkerSize', 8);

hold on;
loglog(iterations_per_wavelength(2:end), relative_error(2:end), 'r+', 'MarkerSize', 10, 'LineWidth', 2.0);
loglog(thit, theory, 'k', 'LineWidth', 2.0);
loglog([0.1, 10E4], relative_error(1)*[1,1], '--','Color',[0,0.4,0.6]);
set(gca,'YMinorGrid','on');
hold off;
legend('Modified Born','PSTD');
fixplot('Iterations per wavelength', 'Relative error', [8 7], '');
xlim([0.1, 1E5]);
ylim([1E-12, 10]);
fp = '../../wavesimpaper/figures/';
print([fp 'accuracy_homogeneous_3d.eps'], '-depsc2');
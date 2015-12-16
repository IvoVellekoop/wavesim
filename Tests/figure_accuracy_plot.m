addpath('..');
load('figure_accuracy.mat');
theory = 1./iterations_per_wavelength.^4 * relative_error(6) * iterations_per_wavelength(6).^4;
loglog(iterations_per_wavelength(2:end), theory(2:end), 'k', 'LineWidth', 2.0);
hold on
loglog([0.1, 10E4], relative_error(1)*[1,1], ':');
loglog(iterations_per_wavelength(1), relative_error(1), 's', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
loglog(iterations_per_wavelength(2:end), relative_error(2:end), 'r+', 'MarkerSize', 10);
xlim([0.1, 5E4]);
hold off;
fixplot('Iterations per wavelength', 'Relative error', [8 7], '');
fp = '../../wavesimpaper/figures/';
print([fp 'accuracy_homogeneous.eps'], '-depsc2');
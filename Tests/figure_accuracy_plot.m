addpath('..');
last = 21;
load('figure_accuracy.mat');
theory = 1./iterations_per_wavelength.^4 * relative_error(2) * iterations_per_wavelength(2).^4;
loglog(iterations_per_wavelength(2:last), theory(2:last), 'k', 'LineWidth', 2.0);
hold on
loglog([0.1, 3E4], relative_error(1)*[1,1], ':');
loglog(iterations_per_wavelength(1), relative_error(1), 's', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
loglog(iterations_per_wavelength(2:last), relative_error(2:last), 'r+', 'MarkerSize', 10);
xlim([0.1, 3E4]);
hold off;
fixplot('Iterations per wavelength', 'Relative error', [8 7], '');
fp = '../../wavesimpaper/figures/';
print([fp 'accuracy_homogeneous.eps'], '-depsc2');
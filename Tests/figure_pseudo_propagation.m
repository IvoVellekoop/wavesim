% Creates the 'pseudo-propagation' figure of the manuscript
addpath('..');

%% options for grid (gopt) and for simulation (sopt) 
PPW=4; %points per wavelength = lambda/h
sopt.lambda = 1; %in mu %lambda_0 = 1; %wavelength in vacuum (in um)
mopt.lambda = sopt.lambda;
mopt.pixel_size = sopt.lambda/PPW;
mopt.boundary_widths = [0, 25*PPW];
mopt.boundary_strength = 0.8;
mopt.boundary_type = 'PML3';

%% define a plane wave source
N = [1, 128*PPW];
sample = SampleMedium(ones(N(1), N(2)), mopt);
source = sparse(N(1), N(2));
source(:,100) = 1; % plane wave source
sim = wavesim(sample, sopt);

fields = zeros(3, N(2));
fields_prev = zeros(3, N(2));
for iteration = 1:3
    it = 20*iteration;
    sim.max_cycles = it/sim.iterations_per_cycle;
    E_wavesim = exec(sim, source);
    fields(iteration, :) = E_wavesim(1,:);
    
    sim.max_cycles = (it-1)/sim.iterations_per_cycle;
    E_wavesim = exec(sim, source);
    fields_prev(iteration, :) = E_wavesim(1,:);
end

fdiff = fields - fields_prev;
xrange = 100:250;
xpos = (xrange-xrange(1))/PPW;
fp = '../../wavesimpaper/figures/';
fdiff = fdiff * 50;
fields = fields * 50;

labels = {'a)', 'b)', 'c)'};
for iteration = 1:3
    it = iteration*20;
    close all;
    figure('Position', [0, 0, 800, 300]);
    plot(xpos, (real(fields(iteration, xrange))));
    %title(it);
    xlim([0, xpos(end)]);
    ylim([-1.1, 1.1]);
    fixplot('x/\lambda', '\psi', [8 3], labels(iteration));
    print([fp 'pseudo_panel_' num2str(iteration) '.eps'], '-depsc2');
end
close all; 
style = {'r', 'b-.', 'k:'};
width = [1.0, 1.0, 1.5];
figure('Position', [0, 0, 800, 400]);
for iteration = 1:3
    it = iteration*20-1;
    plot(xpos, (real(fdiff(iteration, xrange))), style{iteration}, 'LineWidth', width(iteration));
    %title(it);
    xlim([0, xpos(end)]);
    ylimits = [-0.1, 0.1];
    ylim(ylimits);
    set(gca, 'YTick', [-0.1, 0, 0.1]);
    hold on
    plot(it*[1,1]/sim.iterations_per_cycle, ylimits, 'k');
end;
hold off
fixplot('x/\lambda', '\psi', [8 3], 'd)');
print([fp 'pseudo_panel_4.eps'], '-depsc2');
    
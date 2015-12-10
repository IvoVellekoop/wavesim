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

close all;
for iteration = 1:3
    it = iteration*20;
    figure(iteration*2-1);
%    subplot(3, 2, iteration*2-1);
    plot(xpos, (real(fields(iteration, xrange))));
    %title(it);
    xlim([0, xpos(end)]);
    ylim([-0.025, 0.025]);
    xlabel('x/\lambda')
    ylabel('\psi', 'FontSize', 30);
    fixplot();
    save([fp 'pseudo_panel_' num2str(iteration*2-1) '.eps']);
%    subplot(3, 2, iteration*2);
    figure(iteration*2);
    plot(xpos, (real(fdiff(iteration, xrange))), 'r');
    %title(it);
    xlim([0, xpos(end)]);
    ylim([-2E-3, 2E-3]);
    xlabel('x/\lambda');
    ylabel('\psi', 'FontSize', 30);
    hold on
    plot(it*[1,1]/sim.iterations_per_cycle, [-2E-3, 2E-3]);
    hold off
    fixplot();
    save([fp 'pseudo_panel_' num2str(iteration*2) '.eps']);
    drawnow;
end;

% Simulates the propagation of a plane wave through a homogeneous medium
% with n=1
% Analyzes the accuracy of wavesim and PSTD (with varying time steps)

clear all; close all;
WaveSim.use('utilities'); % for analytic solution


%% options for grid and simulation 
PPW=4; %points per wavelength = lambda/h
opt.lambda = 1; %wavelength in vacuum (in um)
opt.callback_interval = 80;
opt.max_iterations = 60000;
opt.energy_threshold = 1E-5;
opt.single_precision = false;

dt_relative_range = 1./2.^(0:0.5:12.5);

opt.pixel_size = opt.lambda/PPW;
opt.boundary_widths = [0, 25*PPW];
opt.boundary_strength = 0.2;
opt.boundary_type = 'PML3'; % PML is deprecated, but we currently don't have anything else for PSTD
N = [1 round(50*PPW)]; % size of medium (in pixels)

%% reserve space for output data
relative_error = zeros(1, length(dt_relative_range)+1);
iterations_per_wavelength = zeros(1, length(dt_relative_range)+1);

%% define a plane wave source and create homogeneous sample
source = Source(ones(N(1),1));
sample = Medium(ones(N), opt);

%% create wavesim object and run the wave propagation simulation
sim = WaveSim(sample, opt);
[E, state] = sim.exec(source);
iterations_per_wavelength(1) = sim.iterations_per_cycle;

%% calculate exact solution analytically
k0 = 2*pi/opt.lambda;
E_theory=homogeneous_medium_analytic_solution(k0, opt.pixel_size, sim.x_range);

% compute relative error of wavesim
difference=E(1,:)-E_theory;
relative_error(1)=mean2(abs(difference).^2) / mean2(abs(E_theory).^2);
disp(relative_error(1))
figure(20);

%% simulate wave propagation for PSTD with varying values for dt
for t_i=1:length(dt_relative_range)   
    % create PSTD object and run the simulation
    opt.dt_relative = dt_relative_range(t_i);
    sim = PSTD(sample, opt);
    iterations_per_wavelength(t_i+1) = sim.iterations_per_cycle;
    [E, state] = sim.exec(source);
    
    %%compare simulation result with exact value
    difference=E(1,:)-E_theory;
    relative_error(t_i+1)=mean2(abs(difference).^2) / mean2(abs(E_theory).^2);
    
    % plot relative errors
    figure(20); clf;
    loglog(iterations_per_wavelength(1), relative_error(1), '+r'); hold on;
    loglog(iterations_per_wavelength(2:t_i+1), relative_error(2:t_i+1), '+b');
    set(gca,'FontSize',14);
    xlabel('Iterations per wavelength','FontSize',14);
    ylabel('Relative error','FontSize',14);
    legend('wavesim','PSTD', 'Location','NorthWest');
end

%save('figure_accuracy.mat','iterations_per_wavelength','relative_error');

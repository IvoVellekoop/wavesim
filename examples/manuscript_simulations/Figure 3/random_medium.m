%%% Simulates the wave propagation of a point source in a random medium
%%% Demonstrate that PSTD converges to wavesim solution as dt decreases

clear all; close all;
addpath('../../../');
addpath('..');

%% options for grid (gopt) and for simulation (opt) 
PPW=4; %points per wavelength = lambda/h
opt.lambda = 1; % wavelength in vacuum (in um)
opt.energy_threshold = 1E-10;
opt.callback_interval = 1000;
opt.max_iterations = 6000;
opt.single_precision = false;

%dt_relative_range = 1/(2^10.5); dt used in manuscript
dt_relative_range = 1/(2^1);

opt.lambda = opt.lambda;
opt.pixel_size = opt.lambda/PPW;
opt.boundary_widths = [0, 0]; %periodic boundaries
opt.boundary_strength = 0;
opt.boundary_type = 'PML3';
N = [64*PPW 64*PPW]; % size of medium (in pixels)

%% Preallocate PSTD data
iterations_per_wavelength = zeros(1, length(dt_relative_range));

%% Construct random medium
% real refractive index
n0 = 1.3;          % mean
n_var = 0.1;     % variance

% imaginary refractive index
a0 = 0.05;       % mean
a_var = 0.02;    % variance

% randomly generate complex refractive index map
n_sample = 1.0*(n0 + n_var * randn(N)) + 1.0i*(a0 + a_var * randn(N));

% low pass filter to remove sharp edges
n_fft = fft2(n_sample);
window = [zeros(1,N(2)/4), ones(1,N(2)/2), zeros(1,N(2)/4)]' * [zeros(1,N(1)/4), ones(1,N(1)/2), zeros(1,N(1)/4)];
n_sample = ifft2(n_fft.*fftshift(window));

% construct sample object
sample = Medium(n_sample, opt); 

%% define a point source at the medium center
source = Source(1, [N(1)/2, N(2)/2]); % point source

%% wavesim simulation
sim = WaveSim(sample, opt);
iterations_per_wavelength(1) = sim.iterations_per_cycle;
E_wavesim = sim.exec(source);

%% PSTD simulation
opt.dt_relative = dt_relative_range;
sim_PSTD = PSTD(sample, opt);
iterations_per_wavelength(2) = sim_PSTD.iterations_per_cycle;
tic;
E_PSTD = sim_PSTD.exec(source);
toc;

error_PSTD = mean2(abs(E_PSTD - E_wavesim).^2) / mean2(abs(E_wavesim).^2);
disp(['Error between PSTD and wavesim: ',num2str(error_PSTD)]);

%% Save results
save('random_medium_results.mat','E_wavesim','N','dt_relative_range','error_PSTD')

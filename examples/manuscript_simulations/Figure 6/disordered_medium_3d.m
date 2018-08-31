%%% Simulates the wave propagation of a point source in a random medium
%%% Demonstrate that PSTD converges to wavesim solution as dt decreases

clear all; close all;
addpath('../../../');
addpath('..');
rng('default'); %reset random number generator

%% options for grid and simulation 
PPW=4; %points per wavelength = lambda/h
opt.lambda = 1; % wavelength in vacuum (in um)
opt.energy_threshold = 0; %limited by max_cycles
opt.callback_interval = 25;
opt.max_cycles = 1100;

dt_relative_range = [0,1/(2^3)];
simulation_run_time = zeros(size(dt_relative_range));

opt.pixel_size = opt.lambda/PPW;
opt.boundary_widths = [0, 0, 0]; %periodic boundaries
N = [32*PPW 32*PPW 32*PPW]; % size of medium (in pixels)

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
n_fft = fftn(n_sample);

W = @(n) [zeros(1,n*3/8), ones(1,n/4), zeros(1,n*3/8)];
window = bsxfun(@times, W(N(2))' * W(N(1)), reshape(W(N(3)), [1,1,N(3)]));
n_sample = ifftn(n_fft.*fftshift(window));
n_sample = max(real(n_sample), 1.0) + 1.0i * max(imag(n_sample), 0.0);

% construct sample object
sample = Medium(n_sample, opt); 

%% define a point source at the medium center
source = Source(1, [N(1)/2, N(2)/2, N(3)/2]); % point source in the center

%% wavesim simulation
sim = WaveSim(sample, opt);
iterations_per_wavelength(1) = sim.iterations_per_cycle;
[E_wavesim, state] = exec(sim, source);
simulation_run_time(1) = state.time;

%% save data
save('disordered_medium_3d.mat','E_wavesim','N','PPW');

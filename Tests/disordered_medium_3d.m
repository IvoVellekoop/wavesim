%%% Simulates the wave propagation of a point source in a random medium
%%% Demonstrate that PSTD converges to wavesim solution as dt decreases

clear all; close all;
addpath('..');
rng('default'); %reset random number generator

%% options for grid (gopt) and for simulation (sopt) 
PPW=4; %points per wavelength = lambda/h
sopt.lambda = 1; %in mu %lambda_0 = 1; %wavelength in vacuum (in um)
sopt.energy_threshold = 1E-30;%16;
sopt.callback_interval = 25;
sopt.max_cycles = 1100;

dt_relative_range = [0,1/(2^3)];
simulation_run_time = zeros(size(dt_relative_range));

mopt.lambda = sopt.lambda;
mopt.pixel_size = sopt.lambda/PPW;
mopt.boundary_widths = [0, 0, 0]; %periodic boundaries
mopt.boundary_strength = 0.2;
mopt.boundary_type = 'PML3';
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
sample = SampleMedium(n_sample, mopt); 

%% define a point source at the medium center
source = zeros(N(1), N(2), N(3));
source(end/2, end/2, end/2) = 1; % point source

%% wavesim simulation
tic;
sim = wavesim(sample, sopt);
iterations_per_wavelength(1) = sim.iterations_per_cycle;
[E_wavesim, state] = exec(sim, source);
simulation_run_time(1) = state.time;

%% PSTD simulations with varying time step size
s2 = sopt;
s2.dt_relative = dt_relative_range(2);
sim_PSTD = PSTD(sample, s2);
iterations_per_wavelength(2) = sim_PSTD.iterations_per_cycle;
[E_PSTD, state] = exec(sim_PSTD, source);
simulation_run_time(2) = state.time;
errors_PSTD = mean(abs(E_PSTD{t_i-1}(:) - E_wavesim(:)).^2) / mean(abs(E_wavesim(:)).^2);

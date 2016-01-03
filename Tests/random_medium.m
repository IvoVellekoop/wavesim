%%% Simulates the wave propagation of a point source in a random medium
%%% Demonstrate that PSTD converges to wavesim solution as dt decreases

clear all; close all;
addpath('..');

%% options for grid (gopt) and for simulation (sopt) 
PPW=4; %points per wavelength = lambda/h
sopt.lambda = 1; %in mu %lambda_0 = 1; %wavelength in vacuum (in um)
sopt.energy_threshold = 1E-25;%16;
sopt.callback_interval = 25;
sopt.max_iterations = 6000;

dt_relative_range = [0,1./2.^(0:0.5:10.5)];

mopt.lambda = sopt.lambda;
mopt.pixel_size = sopt.lambda/PPW;
mopt.boundary_widths = [0, 0]; %periodic boundaries
mopt.boundary_strength = 0;
mopt.boundary_type = 'PML3';
N = [64*PPW 64*PPW]; % size of medium (in pixels)

%% Preallocate PSTD data
E_PSTD = cell(1,length(dt_relative_range)-1);
errors_PSTD = zeros(1,length(E_PSTD));
iterations_per_wavelength = zeros(1, length(dt_relative_range));

%% Construct random medium
% real refractive index
n0 = 1.3;          % mean
n_var = 0.1;     % variance

% imaginary refractive index
a0 = 0.05;        % mean
a_var = 0.02;    % variance

% randomly generate complex refractive index map
n_sample = 1.0*(n0 + n_var * randn(N)) + 1.0i*(a0 + a_var * randn(N));

% low pass filter to remove sharp edges
n_fft = fft2(n_sample);
window = fftshift(hamming(N(1)) * hamming(N(2))');
n_sample = ifft2(n_fft.*window);

% construct sample object
sample = SampleMedium(n_sample, mopt); 

%% define a point source at the medium center
source = sparse(N(1), N(2));
source(end/2,end/2) = 1; % point source

%% wavesim simulation
sim = wavesim(sample, sopt);
iterations_per_wavelength(1) = sim.iterations_per_cycle;
E_wavesim = exec(sim, source);

%% PSTD simulations with varying time step size
for t_i=2:length(dt_relative_range)
    sopt.dt_relative = dt_relative_range(t_i);
    sim_PSTD = PSTD(sample, sopt);
    iterations_per_wavelength(t_i) = sim_PSTD.iterations_per_cycle;
    E_PSTD{t_i-1} = exec(sim_PSTD, source);
    errors_PSTD(t_i-1) = mean2(abs(E_PSTD{t_i-1} - E_wavesim)) / mean2(abs(E_wavesim).^2);
end

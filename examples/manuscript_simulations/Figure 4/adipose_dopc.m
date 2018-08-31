%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 1: point source at focus, used to record forward propagating scattered light
%
% Experiment 2: phase conjugated propagation of result of Experiment 1
% --> forms sharp focus again
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
addpath('../../../');
addpath('..');

%% options for for simulation (opt) and medium (opt)
PPW=4; % points per wavelength
opt.lambda = 1; % wavelength in vacuum (in um)
opt.energy_threshold = 1E-12;
opt.callback_interval = 25;

opt.lambda = opt.lambda;
opt.pixel_size = opt.lambda/PPW;
opt.boundary_widths = [50*PPW, 50*PPW];
opt.boundary_strength = 0.2;
opt.boundary_type = 'PML3';

%% medium properties
%select area of interest
error('Copyrighted image has been removed from the git repository'); 
raw_image = double(imread('adipose_100um_is_52pix_Young2012.png'))/255;
size_per_px = 100/52; %um in image per px

% crop image to square
[a,b] = size(raw_image); s = min(a,b);
raw_image = raw_image(1:s,1:s);

n_min = 1.36; % refractive index corresponding with black pixels
n_max = 1.44; % refractive index corresponding with white pixels
n = raw_image.*(n_max-n_min) + n_min; % linear interpolation based on grayscale image
raw_size = size(raw_image)*size_per_px; %size of refractive index map in um
N = round(raw_size/opt.pixel_size);    %size of refractive index map in pxs
n = imresize(n, N);
sample = Medium(n, opt);

%% Simulations
sim = WaveSim(sample, opt); % wavesim object

%% Experiment 1: point source at desired focus.
% define point source inside scattering medium
tic;
loc_source = [round(N(1) * 3/4), round(N(2)/2 + 1)];
source1 = Source(1, loc_source);

E1 = exec(sim, source1); % perform simulation
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Experiment 2: phase conjugate (top boundary E-field)
% define phase conjugated field at top boundary as source
tic;
source2 = Source(conj(E1(1,:)), [1,1]);

E2 = exec(sim, source2); % perform simulation
toc;

%% Save data
save('dopc_results.mat','E1','E2','opt','n','N');
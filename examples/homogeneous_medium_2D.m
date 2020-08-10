%%% Simulates the wave propagation of a point source in a 2D homogeneous medium
%%% Gerwin Osnabrugge 2015

clear all; close all;
addpath('..');

%% simulation options
PPW=4;                          % points per wavelength
opt.lambda = 1;                 % wavelength in vacuum (in um)
opt.energy_threshold = 1E-10;   % simulation has converged when total added energy is lower than threshold 
opt.pixel_size = opt.lambda/PPW;% grid pixel size (in um)
opt.boundary_widths = PPW*[2,2];% periodic boundaries

%% Construct random medium
N = [32*PPW 32*PPW];     % size of medium (in pixels)
n_sample = ones(N);

%% construct wavesim object
sim = WaveSim(n_sample, opt);

%% define a point source at the center of the medium
source_location = [N(1)/2,N(2)/2];
source = Source(1,source_location);

%% wavesim simulation
[E,state] = exec(sim, source);

%% plot resulting field
figure(1); clf;
imagesc(sim.x_range,sim.y_range,real(E));
axis image;
xlabel('x / \lambda','FontSize',16);
ylabel('y / \lambda','FontSize',16);
h = colorbar;
set(get(h,'Title'),'String','Re\{E\}','FontSize',18,'FontName','Times New Roman');
set(gca,'FontSize',14);
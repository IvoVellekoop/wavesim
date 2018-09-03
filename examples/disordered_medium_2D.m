%%% Simulates the wave propagation of a point source in a 2D random medium
%%% Gerwin Osnabrugge 2015

clear all; close all;
addpath('..');

%% simulations options
PPW=4;                          % points per wavelength
opt.lambda = 1;                 % wavelength in vacuum (in um)
opt.energy_threshold = 1E-8;    % simulation has converged when total added energy is lower than threshold 
opt.pixel_size = opt.lambda/PPW;% grid pixel size (in um)
opt.boundary_widths = [0,0];    % periodic boundaries

%% Construct random medium
% size of medium (in pixels)
N = [64*PPW 64*PPW]; 

% real refractive index
n0 = 1.3;        % mean
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

% construct sample and wavesim object
sample = Medium(n_sample, opt);
sim = WaveSim(sample, opt);

%% define a point source at the center of the medium
source = Source(1,[N(1)/2,N(2)/2]);

%% wavesim simulation
E = exec(sim, source);

%% plot resulting field amplitude
figure(1); clf;

%set axes
x = (-N(2)/2 + 1 : N(2)/2)*opt.pixel_size;
y = (-N(1)/2 + 1 : N(1)/2)*opt.pixel_size;

% plot refractive index distribution
subplot(1,2,1);
imagesc(x,y,real(n_sample));
axis square;
xlabel('x / \lambda','FontSize',16);
ylabel('y / \lambda','FontSize',16);
h = colorbar;
set(get(h,'Title'),'String','n','FontSize',18,'FontName','Times New Roman');
set(gca,'FontSize',14);

% plot resulting field amplitude
subplot(1,2,2);
imagesc(x,y,log(abs(E)));
axis square;
xlabel('x (\lambda)','FontSize',16);
ylabel('y (\lambda)','FontSize',16);
h = colorbar;
set(get(h,'Title'),'String','log|E|','FontSize',18,'FontName','Times New Roman');
set(gca,'FontSize',14);


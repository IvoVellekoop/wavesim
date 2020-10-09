%%% Simulates the wave propagation of a point source in a 3D random medium
%%% Gerwin Osnabrugge 2015

clear all; close all;
addpath('..');

%% simulation options
PPW=4;                           % points per wavelength
opt.lambda = 1;                  % wavelength in vacuum (in um)
opt.energy_threshold = 1E-8;     % simulation has converged when total added energy is lower than threshold 
opt.pixel_size = opt.lambda/PPW; % grid pixel size (in um)
opt.boundary_widths = [0,0,0];   % periodic boundaries

%% Construct Gaussian random medium
% size of medium (in pixels)
N = [32*PPW 32*PPW 32*PPW];     

% real refractive index
n0 = 1.3;        % mean
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

%% construct wavesim object
sim = WaveSim(n_sample, opt);

%% define a point source at the medium center
source = Source(1,[N(1)/2,N(2)/2,N(3)/2]);

%% Perform simulation
E = exec(sim, source);

%% plot resulting field amplitude
figure(1); clf;

% set axes
x = (-N(2)/2 + 1 : N(2)/2)*opt.pixel_size;
y = (-N(1)/2 + 1 : N(1)/2)*opt.pixel_size;
z = (-N(3)/2 + 1 : N(3)/2)*opt.pixel_size;

% plot 3D field amplitude
slice(x,y,z, log(abs(E)), [0 8 16], 17, 17);
set(findobj(gca,'Type','Surface'),'EdgeColor','none')
axis square;
xlabel('x / \lambda','FontSize',16);
ylabel('y / \lambda','FontSize',16);
zlabel('z / \lambda','FontSize',16);
h = colorbar;
set(get(h,'Title'),'String','log|E|','FontSize',18,'FontName','Times New Roman');
set(gca,'FontSize',14);
xlim([0 16]);


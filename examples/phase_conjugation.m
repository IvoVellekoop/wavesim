%%% simulates phase conjugation experiment in random medium
%%% Experiment 1: point source at focus, used to record forward propagating scattered light
%
%%% Experiment 2: phase conjugated propagation of result of Experiment 1
% --> forms sharp focus again
%%% Gerwin Osnabrugge 2015

clear all; close all;
addpath('..');

%% simulations options
PPW=4;                                    % points per wavelength
opt.lambda = 1;                           % wavelength in vacuum (in um)
opt.energy_threshold = 1E-12;             % simulation has converged when total added energy is lower than threshold 
opt.pixel_size = opt.lambda/PPW;          % grid pixel size (in um)
opt.boundary_widths = [4*PPW,4*PPW];      % absorbing boundaries
opt.usemex = true;
if(opt.usemex)
    addpath('..\MexBin');
end
%% Construct random medium
% size of medium (in pixels)
N = PPW*[128 64];

% real refractive index
n0 = 1.3;          % mean
n_var = 0.1;     % variance

% randomly generate complex refractive index map
n_sample = 1.0*(n0 + n_var * randn(N));

% low pass filter to remove sharp edges
n_fft = fft2(n_sample);
cutoff = 1/16; % 1 is no cutoff frequency
window = ([zeros(1,N(2)*(1/2-cutoff)), ones(1,N(2)*(2*cutoff)), zeros(1,N(2)*(1/2-cutoff))]' * ...
    	 [zeros(1,N(1)*(1/2-cutoff)), ones(1,N(1)*(2*cutoff)), zeros(1,N(1)*(1/2-cutoff))])';
n_sample = real(ifft2(n_fft.*fftshift(window)));

%% construct wavesim object
sim = WaveSim(n_sample, opt);

%% Experiment 1: recording phase
% point source in random medium
source_position = [N(1)*3/4,N(2)*1/2];
source1 = Source(1,source_position);    

% run simulation
E_recording = exec(sim, source1);

%% Experiment 2: playback phase
% use phase conjugated field at input plane as source
Econj = conj(E_recording(1,:));     % conjugated field at input plane
source2 = Source(Econj,[1,1]);

% run simulation
E_playback = exec(sim, source2);

%% plot resulting field amplitude
fig = figure(1); clf;
set(fig,'Position',get(fig,'Position')+[0,0, 400, 0]);

% set axes
x = (-N(1)/2 + 1 : N(1)/2)*opt.pixel_size;
y = (-N(2)/2 + 1 : N(2)/2)*opt.pixel_size;

% plot recording phase
subplot(1,2,1);
imagesc(x,y,log(abs(E_recording)));
title('Recording Phase','FontSize',16);
axis square;
xlabel('x / \lambda','FontSize',16);
ylabel('y / \lambda','FontSize',16);
h = colorbar;
set(get(h,'Title'),'String','log|E|','FontSize',18,'FontName','Times New Roman');
set(gca,'FontSize',14);

% plot playback phase
subplot(1,2,2);
imagesc(x,y,log(abs(E_playback)));
axis square;
title('Playback Phase','FontSize',16);
xlabel('x / \lambda','FontSize',16);
ylabel('y / \lambda','FontSize',16);
h = colorbar;
set(get(h,'Title'),'String','log|E|','FontSize',18,'FontName','Times New Roman');
set(gca,'FontSize',14);


%%% Simulate the pulse propagation of a plane wave with a specified bandwidth 
%%% in a 2D disordered medium

clear all; close all;
addpath('..');

%% simulations options
PPW=4;                            % points per wavelength
opt.lambda = 1;                   % wavelength in vacuum (in um)
opt.energy_threshold = 1E-10;     % simulation has converged when total added energy is lower than threshold 
opt.pixel_size = opt.lambda/PPW;  % grid pixel size (in um)
opt.boundary_widths = [0,4*PPW];  % periodic boundary in y, absorbing boundary in x

%% Construct homogeneous medium
N = [32*PPW 32*PPW];            % size of medium (in pixels)
n_sample = ones(N);             % sample refractive index

%% define a source with spectrum (800-1200nm) at the medium center
lambda = 1;                                 % center wavelength (in um)
Nfreq = 50;                                 % number of independent frequencies simulated
lambda_set = linspace(0.8,1.2,Nfreq)*lambda;% source bandwidth
source_spectrum = gausswin(Nfreq,4);          % source spectrum

%% wavesim simulation for all source frequencies
% Preallocate data for fields
E_set = zeros(N(1), N(2), Nfreq);

for f = 1:Nfreq
    % change source wavelength and amplitude
    disp(['Simulation: ',num2str(f),'/',num2str(Nfreq)]);
    source = Source(ones(N(1),1)*source_spectrum(f),[1,1]); %plane wave
    opt.lambda = lambda_set(f);
    
    % run simulation and store resulting field
    sim = WaveSim(n_sample, opt);
    E_set(:,:,f) = exec(sim, source);
end

%% field in time domain
% Inverse Fourier transform E(frequency) -> E(time)
E_time = ifftshift(ifft(E_set,Nfreq,3),3);

% set axes
x = (1 : N(2))*opt.pixel_size;
y = (-N(1)/2 + 1 : N(1)/2)*opt.pixel_size;
f = 3e8*2*pi./(lambda_set*10^-6);                  % frequency: f = ck0 [in Hz]
t_set = (-Nfreq/2+1:Nfreq/2)* 1/( 2*(f(1)-f(2)) );

% only store signal where time is positive
E_time = E_time(:,:,t_set >= 0);
t_set = t_set(t_set >= 0) * 10^15;       % in fs

% plot all fields as funcion of time
figure(2); clf;

Imax = max(abs(E_time(:)).^2);
for t = 1:Nfreq/2
    figure(2); imagesc(x,y,abs(E_time(:,:,t)).^2,[0,Imax]);
    title(['t = ',num2str(t_set(t)),' fs']);
    xlabel('x / \lambda','FontSize',16);
    ylabel('y / \lambda','FontSize',16);
    h = colorbar;
    set(get(h,'Title'),'String','log|E|','FontSize',18,'FontName','Times New Roman');
    set(gca,'FontSize',16);
    pause(0.1);    
end
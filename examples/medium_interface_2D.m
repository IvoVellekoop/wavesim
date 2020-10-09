%%% example script simulating a plane wave through an interface of two
%%% media with different refractive indices

clear
%% simulation options
PPW=8;                           % points per wavelength
opt.lambda = 1;                  % wavelength in vacuum (in um)
opt.energy_threshold = 1E-4;     % simulation has converged when total added energy is lower than threshold 
opt.pixel_size = opt.lambda/PPW; % grid pixel size (in um)
opt.boundary_widths = PPW*[2,2]; % absorbing boundaries

%% sample properties
N = PPW*[16,32]; % size of medium (in pixels)
n1 = 1;          % refractive index medium 1 (top)
n2 = 2;          % refractive index medium 2 (bottom)

%% make refractive index distribution and create WaveSim object
n_sample = [n1*ones(N(1)/2,N(2));n2*ones(N(1)/2,N(2))];
sim = WaveSimVector(n_sample, opt);

%% define plane wave source with Gaussian intensity profile with incident angle theta
% properties
theta = pi/4;                            % angle of plane wave
kx = 2*pi/opt.lambda*sin(theta);
x = (1:N(2))*opt.pixel_size;             

% create source object
source_amplitude = [gausswin(N(2)/2,3)'.*exp(1.0i*kx*x(1:end/2)),zeros(1,N(2)/2)];
source = Source(source_amplitude, [1,1,1,1]); % horizontally-polarized plane wave source under a angle

%% perform simulation and plot resulting field
E = sim.exec(source);

% plot results (in horizontal-direction only)
figure;
imagesc(sim.x_range,sim.y_range,abs(E(:,:,:,1)));
xlabel('x / \lambda','FontSize',16);
ylabel('y / \lambda','FontSize',16);
h = colorbar;
set(get(h,'Title'),'String','|E|','FontSize',18,'FontName','Times New Roman');
set(gca,'FontSize',14);
axis image;

%% Simulates light scattering on a Mie scatterer

function Escat = run_mie_sphere_simulation(opt, n_sphere, z_r)
    %% create wavesim object
    sim = WaveSimVector(n_sphere, opt);

    %% Define source.    
    % calculate source prefactor (see analytical solution for homogeneous medium)
    k = opt.bg_index*2*pi/opt.lambda;
    prefactor = 1.0i*opt.pixel_size/(2*k);
    
    % source intensity profile (filter edges to reduce diffraction and the
    % source extends into the absorbing boundaries)
    sx = size(n_sphere, 1);
    sy = size(n_sphere, 2);
    srcx = reshape(tukeywin(sx, 0.5), [1, sx, 1]);
    srcy = reshape(tukeywin(sy, 0.5), [sy, 1, 1]);
    source_amplitude = 1/prefactor * exp(1.0i*k*z_r(1)) * srcx .* srcy;
    
    % create source object  
    source = Source(source_amplitude, [1, 1, 1, 1]);
    
    %% run Mie sphere simulation
    E = sim.exec(source);

    %% run similar simulation, now without medium (to get background wave)
    opt.epsilon = sim.epsilon;
    
    n_sample = opt.bg_index*ones(size(n_sphere));
    sim = WaveSimVector(n_sample, opt);
    
    Ebg = sim.exec(source);
    Escat = E-Ebg;    
end
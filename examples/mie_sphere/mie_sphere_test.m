%%% script used to simulate a mie sphere (non trivial solution) with 
%%% different boundary conditions. Simulation results are compared with
%%% solution obtained from Mie theory

clear; close all; clc; warning off;
addpath('../../'); % add path to wavesim core code

%% load analytical result for Mie sphere
try 
    load('mie_sphere_theory.mat');
catch
    error('Unable to read file: mie_sphere_theory.mat. First run mie_sphere_analytical_solution to obtain analytical results');
end

%% simulations settings
opt.energy_threshold    = 1E-99;            % no energy threshold to make sure the same of iterations are performed in all simulations
opt.max_cycles          = 200;              % maximum number of wave periods to run the simulation
opt.single_precision    = true;

% callback settings
opt.callback_interval   = 8;
opt.callback            = @Simulation.abs_crossimage_callback;

% tested parameters
boundary_widths = (1:12)';                  % tested boundary widths (in wavelengths)
boundary_labels = {'PBL','ARL','PBL (with ACC)','ARL (with ACC)'}; % boundary type

%% Simulate Mie scatternig for window and PML3 boundaries with varying boundary widths
% initialize result dataset
sim_results = zeros(numel(boundary_widths), numel(boundary_labels));
E_set = cell(numel(boundary_widths), numel(boundary_labels));

% calculate simulation results
for w = 1:numel(boundary_widths)
    % change boundary width
    opt.boundary_widths = boundary_widths(w)*ones(1,3)/opt.pixel_size; 
    for b = 1:4
        % change boundary type
        opt = boundary_type(opt,b);
        
        % run simulation (including background simulation)
        E_sim = run_mie_sphere_simulation(opt, n_sphere, z_r);
        
        % calculate relative error
        sim_results(w,b) = mean(abs(E_sim(:) - E_theory(:)).^2) / mean(abs(E_theory(:)).^2);
        
        % store suimulated field
        E_set{w,b} = E_sim;
    end
end

%% separately store field cross sections for manuscript figure
% xz-cross section at y=0 with boundary width of 6 lambda

% x-polarization
Ex_theory = squeeze(E_theory(end/2,:,:,1));
Ex_PBL = squeeze(E_set{6,1}(end/2,:,:,1));
Ex_ARL = squeeze(E_set{6,4}(end/2,:,:,1));

% z-polarization
Ez_theory = squeeze(E_theory(end/2,:,:,3));
Ez_PBL = squeeze(E_set{6,1}(end/2,:,:,3));
Ez_ARL = squeeze(E_set{6,4}(end/2,:,:,3));

%% save simulation results
% save('mie_sphere_results.mat','boundary_labels','boundary_widths',...
%      'Ex_ARL','Ex_PBL','Ex_theory','Ez_ARL','Ez_PBL','Ez_theory',...
%      'sim_results','x_r','y_r','z_r');

%% function with all tested types of boundary conditions
function opt = boundary_type(opt,boundary_number)
switch boundary_number
    case 1 % PBL3 boundary without ACC
        opt.boundary_type = 'PBL3';
        opt.boundary_strength = 0.2;
        opt.ACC = false;
    case 2 % window boundary without ACC
        opt.boundary_type = 'ARL';
        opt.ACC = false;
    case 3 % PBL3 boundary with ACC
        opt.boundary_type = 'PBL3';
        opt.boundary_strength = 0.2;
        opt.ACC= true;
    case 4 % window boundary with ACC
        opt.boundary_type = 'ARL';
        opt.ACC = true;
    otherwise
        error(['boundary number ',num2str(boundary_number),' is not defined']);
end
end
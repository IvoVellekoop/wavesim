classdef Simulation
    % Base class for a 2-D or 3D wave simulation
    % Ivo Vellekoop & Gerwin Osnabrugge 2016-2020
    %
    % Data storage:
    % Internally, all data is stored as 2-D or 3-D arrays. For vector simulations,
    % the data for each polarization components is stored as a separate
    % array in a 1x3 cell array.
    %
    % The size of the simulation grid is determined by the refractive index
    % distribution that is passed in the constructor (see Medium.m). 
    % The grid size includes the region of interest (ROI), boundaries, and
    % padding for efficient application of the fft algorithm. In the end,
    % only the field within the ROI is returned by exec(), although it is
    % possible to recover the full field including boundaries and padding
    % through the 'state' variable that is an optional return argument.
    %
    properties
        sample;   % a medium object containing the map of the relative 
        % dielectric constant including the boundary layers
        %options:
        roi = []; % Part of the simulated field that is returned as output,
        % this matrix contains 2 rows, one for the start index and one for
        % the last index. It contains 4 columns (one for each dimension and
        % one to indicate the polarization component, which is always 1 for
        % scalar simulations)
        % Defaults to the full medium. To save memory, a smaller
        % output_roi can be specified.         
        lambda = 1;      % Wavelength (in micrometers)        
        gpu_enabled = gpuDeviceCount > 0; % flag to determine if simulation are run
        % on the GPU (default: run on GPU if we have one)        
        single_precision = true; % flag to determine if single precision or
        % double precision calculations are used.
        % Note that on a typical GPU, double
        % precision calculations are about 10
        % times as slow as single precision.
        
        % Callback function options
        % The simulation code calls a callback function every few
        % iteratations. There are two callbacks already implemented:
        % default_callback and abs_image_callback.
        % you can choose which one to use to set the 'callback' option
        % note that you can also specify your own custom callback function
        callback = @Simulation.default_callback; % callback function that is called for showing the progress of the simulation. Default shows image of the absolute value of the field.
        callback_options = struct; %other options for the callback
        callback_interval = 50; % the callback is called every 'callback_interval' iterations.
        
        % Stopping criteria
        % To determine when to stop, the algorithm calculates the amount of
        % energy added in the last iteration. The
        % simulation is stopped when the total energy in this 'Ediff' is
        % lower than a threshold value.
        %
        energy_threshold = 1E-2; % Threshold for terminating the simulation.
        % The threshold is specified as a fraction
        % of the amount of energy added during the
        % first iteration of the algorithm. Therefore
        % it is compensates for the
        % amplitude of the source terms.
        
        energy_calculation_interval = 8; % only calculate energy difference every N steps (to reduce overhead)
        max_cycles = inf; % Maximum number of wave periods to run the simulation.
        % Note that the number of actual iterations per optical cycle
        % depends on the algorithm and its parameters
        % (notably the refractive index contrast).
        % Therefore, it is not recommended to specify
        % max_cycles explicitly.
        
        %internal:
        grid; % simgrid object
        N; % size of simulation grid (in pixels)
        x_range; %x-coordinates of roi
        y_range; %y-coordinates of roi
        z_range; %z-coordinates of roi
        %internal:
        iterations_per_cycle;%must be set by derived class
    end
    
    methods
        function obj = Simulation(refractive_index, options)
            %% Constructs a simulation object
            %	sample = SampleMedium object
            %   options = simulation options.
            %   options.max_cycles = maximum number of optical cycles for
            %   which to run the simulation. Note that the number of
            %   required iterations per cycle depends on the algorithm and
            %   its parameters
     
            % copy values from 'options' to 'obj' (only if the field is present in obj)
            fs = intersect(fields(obj), fields(options));
            for f=1:length(fs)
                fieldname = fs{f};
                obj.(fieldname) = options.(fieldname);
            end
            
            % generate Medium object and set corresponding grid
            [obj.sample, obj.grid] = Medium(refractive_index, options);
            
            % set grid and number of grid points
            obj.N    = [obj.grid.N, 1]; %vector simulations set last parameter to 3
            
            % set region of interest
            % by default the full field is returned (excluding boundaries)
            default_roi = [obj.sample.Bl + 1; obj.grid.N - obj.sample.Br];
            if isempty(obj.roi)     
                obj.roi = [Simulation.make4(default_roi(1,:)); Simulation.make4(default_roi(2,:))]; 
            else % use specified region of interest
                obj.roi = [Simulation.make4(obj.roi(1,:)); Simulation.make4(obj.roi(2,:))];
                obj.roi(:,1:3) = obj.roi(:,1:3) + default_roi(1,:) - 1; % shift xyz to match grid coordinates
            end
            
            % calculate roi coordinates
            obj.x_range = obj.grid.x_range(obj.roi(1,2):obj.roi(2,2));
            obj.x_range = obj.x_range - obj.x_range(1);
            obj.y_range = obj.grid.y_range(obj.roi(1,1):obj.roi(2,1));
            obj.y_range = obj.y_range - obj.y_range(1);
            obj.z_range = obj.grid.z_range(obj.roi(1,3):obj.roi(2,3));
            obj.z_range = obj.z_range - obj.z_range(1);
            
            %update the energy at least every callback
            obj.energy_calculation_interval = min(obj.energy_calculation_interval, obj.callback_interval); 
        end
        
        function [E, state] = exec(obj, source)
            % EXEC Executes the simulation with the given source
            % source    is a 'source' object (see documentation for source)
            %
            tic;
            % convert the source to the correct data type (single, double,
            % gpuarray or not).
            % shift source so that [1,1,1] corresponds to start of roi
            % then crops source so that anything outside the simulation
            % grid is removed. Note: it is allowed to put a
            % source energy outside of the roi, but inside the grid (that
            % is, inside the absorbing boundary).
            %
            if (~isa(source, 'Source'))
                error('Expecting a Source object as input parameter');
            end
            state.source = obj.data_array(source)...
                              .shift(obj.roi(1,:))...
                              .crop([1,1,1,1; obj.N]);
            state.source_energy = sum(state.source.energy);
            if state.source_energy == 0
                warning('There is no source inside the grid boundaries, expect to get zero field everywhere.');
                return;
            end
            
            %%% prepare state (contains all data that is unique to a single
            % run of the simulation)
            %
            state.it = 1; %iteration
            state.max_iterations = ceil(obj.max_cycles * obj.iterations_per_cycle);
            state.diff_energy = [];
            state.last_step_energy = inf;
            state.calculate_energy = true;
            state.has_next = true;
            
            %%% Execute simulation
            % the run_algorithm function should:
            % - update state.last_step_energy when required
            % - call next(obj, state)
            % - return calculated field in state.E
            %(see wavesim for an example)
            state = run_algorithm(obj, state);
            
            state.time=toc;
            if state.converged
                disp(['Reached steady state in ' num2str(state.it) ' iterations']);
            else
                disp('Did not reach steady state');
            end
            disp(['Time consumption: ' num2str(state.time) ' s']);
            
            E = gather(state.E); % convert gpuArray back to normal array
        end
        
        %% Continue to the next iteration. Returns false to indicate that the simulation has terminated
        function state = next(obj, state, can_terminate)
            %% store energy (relative to total energy in source)
            state.diff_energy(state.it) = state.last_step_energy;
            
            %% check if simulation should terminate
            if can_terminate && state.diff_energy(state.it) / state.diff_energy(1) < obj.energy_threshold
                state.has_next = false;
                state.converged = true;
            elseif can_terminate && state.it >= state.max_iterations
                state.has_next = false;
                state.converged = false;
            else
                state.has_next = true;
            end
            
            %% do we need to calculate the last added energy? (only do this once in a while, because of the overhead)
            state.calculate_energy = mod(state.it - 1, obj.energy_calculation_interval)==0;
            
            %% call callback function if neened
            if (mod(state.it, obj.callback_interval)==0 || ~state.has_next) %now and then, call the callback function to give user feedback
                obj.callback(obj, state);
            end
            state.it = state.it+1;
        end
    end
    
    methods(Access = protected)
        function d = data_array(obj, data, N)
            % data_array(obj, data) - converts the data to an array of the correct format
            %                         if gpu_enabled is true, the array is created on the gpu
            %                         the data is converted to
            %                         single/double precision based on the
            %                         options set when creating the
            %                         simulation object.
            % data_array(obj, [], N)- creates an empty array (filled with zeros) of size N
            % Check whether single precision and gpu computation options are enabled
            if obj.single_precision
                p = 'single';
            else
                p = 'double';
            end
            
            %obj.N
            if isempty(data) %no data is specified, generate empty array of specified size
                if obj.gpu_enabled
                    d = zeros(N, p, 'gpuArray');
                else
                    d = zeros(N, p);
                end
                return;
            end
            d = cast(full(data), p);
            if obj.gpu_enabled
                d = gpuArray(d);
            end
        end
        
        function Ecrop = crop_field(obj,E)
            % Removes the boundary layer from the simulated field by
            % cropping field dataset
            Ecrop = E(obj.roi(1,1):obj.roi(2,1),...
                      obj.roi(1,2):obj.roi(2,2),...
                      obj.roi(1,3):obj.roi(2,3),...
                      obj.roi(1,4):obj.roi(2,4));
        end
    end
    
    methods(Static) 
        function en = energy(E_x, roi)
            % calculculate energy (absolute value squared) of E_x
            % optionally specify roi to indicate which values are to be
            % included
            if nargin == 2
                E_x = E_x(roi(1,1):roi(2,1), roi(1,2):roi(2,2), roi(1,3):roi(2,3), roi(1,4):roi(2,4));
            end
            en = full(gather(sum(abs(E_x(:)).^2)));
        end
        
       % callback functions
        function default_callback(obj, state, dim, pol)
            %default callback function. Shows  total energy evolution and 
            %real value of field along specified dimension (longest
            %dimension by default)
            %
            figure(1); clf;
            energy = state.diff_energy(1:state.it) / state.diff_energy(1);
            threshold = obj.energy_threshold;
            E = state.E;
            subplot(2,1,1); plot(1:length(energy),log10(energy),'b',[1,length(energy)],log10(threshold)*ones(1,2),'--r');
            title(length(energy));  xlabel('# iterations'); ylabel('log_{10}(energy added)');
            
            % plot midline along longest dimension
            xpos = ceil(size(E, 2)/2);
            ypos = ceil(size(E, 1)/2);
            zpos = ceil(size(E, 3)/2);
            if nargin < 3
                [~,dim] = max([ypos,xpos,zpos]);
            end
            if nargin < 4 % default show horizontal polarization
                pol = 1;
            end
            
            switch dim
                case 1 % y-dimension
                    sig = log(abs(E(:, xpos, zpos,pol)));
                    ax = obj.y_range(:);
                    sroi = ax(obj.roi(:,1));
                    ax_label = 'y (\mum)';
                case 2 % x-dimension
                    sig = log(abs(E(ypos, :, zpos,pol)));                    
                    ax = obj.x_range(:);
                    sroi = ax(obj.roi(:,2));
                    ax_label = 'x (\mum)';
                case 3 % z-dimension
                    sig = log(abs(E(ypos, xpos, :,pol)));   
                    ax = obj.grid.z_range(:);
                    sroi = ax(obj.roi(:,3));
                    ax_label = 'z (\mum)';
            end
            subplot(2,1,2); 
            plot(ax, squeeze(sig),'b'); hold on;
            title('midline cross-section')
            xlabel(ax_label); ylabel('log(|E|)');
            
            % draw dashed lines to indicate start of boundary layer
            ca = gca;
            plot(sroi(1)*ones(1,2), [ca.YLim(1),ca.YLim(2)],'--k', ...
                 sroi(end)*ones(1,2), [ca.YLim(1),ca.YLim(2)],'--k');
            
            %disp(['Added energy ', num2str(energy(end))]);
            drawnow;
        end
        
        function abs_crossimage_callback(obj, state)
            % callback function that displays the intensity along the
            % propagation direction
            %
            figure(1);
            if size(state.E,3) == 1 % 1D/2D simulation
                Eprop = abs(squeeze(state.E(:,:,1,1)));
            else
                Eprop = abs(squeeze(state.E(ceil(end/2),:,:,1)));
            end
            imagesc(Eprop);
            axis image;
            title(['Differential energy ' num2str(state.diff_energy(state.it)/state.diff_energy(1))]);
            drawnow;
        end
        
        %skip the callback function, for benchmarking speed
        function no_callback(obj, state)
        end
        
        function sz = make3(sz, pad)
            % makes sure sz is a row vector with exactly 3 elements
            % (pads when needed)
            % (this function is needed because 'size' by default removes
            %  trailing singleton dimensions)
            sz = sz(:).'; %make row vector
            if numel(sz) < 3
                sz((end+1):3) = pad;
            end
        end
        
        function sz = make4(sz)
            % converts the vector sz to a 4-element vector by appending 1's
            % when needed. This function is useful since 'size' removes
            % trailing singleton dimensions of arrays, so a 100x100x1x1
            % array returns a size of [100, 100], whereas a 100x100x1x3
            % array retuns a size of [100, 100, 1, 3]
            % As a workaround for this inconsistency, we always use
            % 4-element size vectors.
            sz = sz(:).';
            if numel(sz) < 4
                sz((end+1):4) = 1;
            else
                if numel(sz) > 4
                    error('Position and size vectors can be 4-dimensional at most');
                end
            end
        end
        
        % e.g. call Simulation.use('utilities') to add the 'utilities' folder
        % to the MATLAB path
        function use(subpath)
            addpath([fileparts(mfilename('fullpath')) '\' subpath]); %todo: put in more sensible location
        end
    end
end
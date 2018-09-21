classdef Simulation
    %Base class for a 2-D wave simulation
    % Ivo M. Vellekoop 2015-2018
    %
    % Data storage:
    % Internally, all data is stored as 2-D or 3-D arrays. For vector simulations,
    % the data for each polarization components is stored as a separate
    % array in a 1x3 cell array.
    %
    % The size of the simulation grid is determined by the sample object
    % that is passed in the constructor (see SampleMedium.m). The grid size
    % includes the region of interest (ROI), boundaries, and
    % padding for efficient application of the fft algorithm. In the end,
    % only the field within the ROI is returned by exec(), although it is
    % possible to recover the full field including boundaries and padding
    % through the 'state' variable that is an optional return argument.
    %
    properties
        %options:
        output_roi = []; % Part of the simulated field that is returned as output, 
                         % Defaults to the full medium. To save memory, a smaller
                         % output_roi can be specified. [experimental]
        
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
        
        % Stopping options
        % To determine when to stop, the algorithm calculates the
        % difference in the field before and after an iteration. The
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
        roi; % position of simulation area with respect to padded array.
        %      this matrix contains 2 rows, one for the start index and one for
        %      the last index. It contains 4 columns (one for each dimension and
        %      one to indicate the polarization component, which is always 1 for
        %      scalar simulations)
        N; % size of simulation grid (in pixels)
        x_range; %x-coordinates of roi
        y_range; %y-coordinates of roi
        z_range; %z-coordinates of roi
        %internal:
        iterations_per_cycle;%must be set by derived class
    end
    
    methods
        function obj = Simulation(sample, options)
            %% Constructs a simulation object
            %	sample = SampleMedium object
            %   options = simulation options. 
            %   options.max_cycles = maximum number of optical cycles for
            %   which to run the simulation. Note that the number of
            %   required iterations per cycle depends on the algorithm and
            %   its parameters

            %copy values from 'options' to 'obj' (only if the field is present in obj)
            fs = intersect(fields(obj), fields(options));
            for f=1:length(fs)
                fieldname = fs{f};
                obj.(fieldname) = options.(fieldname);
            end
            
            obj.grid = sample.grid;
            obj.roi  = [sample.roi, [1;1]]; %vector simulation objects change the last column to [1;3] to indicate 3 polarization channels.
            if isempty(obj.output_roi)
                obj.output_roi = obj.roi; % by default return the full field
            else
                obj.output_roi = [Source.make4(obj.output_roi(1,:)); Source.make4(obj.output_roi(2,:))];
                obj.output_roi = obj.output_roi + obj.roi(1,:) - 1; %shift to match grid coordinates
            end
            
            obj.energy_calculation_interval = min(obj.energy_calculation_interval, obj.callback_interval); %update the energy at least every callback
                
            obj.N    = [obj.grid.N, 1]; %vector simulations set last parameter to 3
            obj.x_range = sample.grid.x_range(obj.roi(1,2):obj.roi(2,2));
            obj.x_range = obj.x_range - obj.x_range(1);
            obj.y_range = sample.grid.y_range(obj.roi(1,1):obj.roi(2,1));
            obj.y_range = obj.y_range - obj.y_range(1);
            obj.z_range = sample.grid.z_range(obj.roi(1,3):obj.roi(2,3));
            obj.z_range = obj.z_range - obj.z_range(1);
            
            % determine how many optical cycles to simulate.
            % By default, simulate enough cycles to pass the medium 1.5
            % times. Note that this will not be sufficient if there are
            % significant reflections.
            %
            %if (obj.max_cycles == 0)
            %    LN = sqrt(length(obj.x_range).^2 + length(obj.y_range)^2 + length(obj.z_range)^2);
            %    obj.max_cycles = LN * obj.grid.dx / obj.lambda * 2;
            %end
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
                disp(['Time consumption: ' num2str(state.time) ' s']);
            else
                disp('Did not reach steady state');
                disp(['Time consumption: ' num2str(state.time) ' s']);
            end
            
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
        
        function abs_image_callback(obj, state)
            % by default this function displays state.E, but it can be
            % configured to show different variables, such as state.Ediff
            %
            figure(1);
            if isfield(obj.callback_options, 'fieldname')
                E = state.(obj.callback_options.fieldname);
            else
                E = state.E;
            end
            imagesc(abs(E(:,:,ceil(end/2), 1, ceil(end/2))));
            axis image;
            title(['Differential energy ' num2str(state.diff_energy(state.it)/state.diff_energy(1))]);
            drawnow;
        end
        
        %skip the callback function, for benchmarking speed
        function no_callback(obj, state)
        end
        
        %default callback function. Shows real value of field, and total energy evolution
        function default_callback(obj, state)
            figure(1);
            energy = state.diff_energy(1:state.it) / state.diff_energy(1);
            threshold = obj.energy_threshold;
            E = state.E;
            subplot(2,1,1); plot(1:length(energy),log10(energy),'b',[1,length(energy)],log10(threshold)*ones(1,2),'--r');
            title(length(energy));  xlabel('# iterations'); ylabel('log_{10}(energy added)');
            
            % plot midline along longest dimension
            xpos = ceil(size(E, 2)/2);
            ypos = ceil(size(E, 1)/2);
            zpos = ceil(size(E, 3)/2);
            if xpos > max(ypos, zpos)
                sig = log(abs(E(ypos, :, zpos)));
                sroi = obj.roi(:,1);
            elseif ypos > max(xpos, zpos)
                sig = log(abs(E(:, xpos, zpos)));
                sroi = obj.roi(:,2);
            else
                sig = log(abs(E(ypos, xpos, :)));
                sroi = obj.roi(:,3);
            end
            subplot(2,1,2);
            plot(1:length(sig), squeeze(sig), sroi(1)*ones(1,2), [min(sig),max(sig)], sroi(end)*ones(1,2), [min(sig),max(sig)]);
            title('midline cross-section')
            xlabel('y (\lambda / 4)'); ylabel('log(|E_x|)');
            
            %disp(['Added energy ', num2str(energy(end))]);
            drawnow;
        end
    end
end
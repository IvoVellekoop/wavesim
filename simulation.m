classdef simulation
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
        lambda = 1; %wavelength (default = 1 unit)      
        gpu_enabled = gpuDeviceCount > 0; % flag to determine if simulation are run on the GPU (default: run on GPU if we have one)
        single_precision = false; % flag to determine if single precision is used (default: double)
        callback = @simulation.default_callback; %callback function that is called for showing the progress of the simulation. Default shows image of the absolute value of the field.
        callback_interval = 50; %the callback is called every 'callback_interval' steps. Default = 5
        energy_threshold = 1E-20; %the simulation stops when the difference for a step is less than 'energy_threshold'
        energy_calculation_interval = 10; %only calculate energy difference every N steps (to reduce overhead)
        max_cycles = 0; %number of wave periods to run the simulation. The number of actual iterations per cycle depends on the algorithm and its parameters
        %internal:
        iterations_per_cycle;%must be set by derived class
    end
    
    methods
        function obj=simulation(sample, options)
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
            obj.N    = [obj.grid.N, 1]; %vector simulations set last parameter to 3
            obj.x_range = sample.grid.x_range(obj.roi(1,1):obj.roi(2,1));
            obj.x_range = obj.x_range - obj.x_range(1);
            obj.y_range = sample.grid.y_range(obj.roi(1,2):obj.roi(2,2));
            obj.y_range = obj.y_range - obj.y_range(1);
            obj.z_range = sample.grid.z_range(obj.roi(1,3):obj.roi(2,3));
            obj.z_range = obj.z_range - obj.z_range(1);
            
            % determine how many optical cycles to simulate.
            % By default, simulate enough cycles to pass the medium 1.5
            % times. Note that this will not be sufficient if there are
            % significant reflections.
            %
            if (obj.max_cycles == 0)
                LN = sqrt(length(obj.x_range).^2 + length(obj.y_range)^2 + length(obj.z_range)^2);
                obj.max_cycles = LN * obj.grid.dx / obj.lambda * 2;
            end
        end
        
        function [E, state] = exec(obj, source, source_pos)
 %% exec Executes the simulation with the given source
% source    is a 2-D or 3-D array with source data.
%           The size of 'source' should not exceed the size of
%           the roi of the simulation grid, but it can be
%           smaller. In that case the coordinate vector 
%           source_pos can be used to indicate the starting
%           index of the source within the refractive index map
% source_pos  (optional for scalar simulations, defaults to [1,1,1]), can be specified to
%           position the source with respect to the top left corner of the
%           simulation grid. For vector simulations, this vector contains
%           one extra coordinate (e. g. [1,1,1,1] for a 3-D simulation)
%           to indicate the polarization of the source.
%
% It is possible to specify multiple sources by passing cell arrays for
% 'sources' and 'positions'
%
            tic;
            if nargin < 3
                source_pos = [1,1,1];
            end
            [state.source, state.source_pos, state.source_energy] = prepare_source(obj, source, source_pos);
            
            %%% prepare state (contains all data that is unique to a single 
            % run of the simulation)
            %
            state.it = 1; %iteration
            state.max_iterations = ceil(obj.max_cycles * obj.iterations_per_cycle);
            state.diff_energy = zeros(state.max_iterations,1);
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
            
            %% return only part inside roi. Array remains on the gpu if gpuEnabled = true
            E = state.E(obj.roi(1,1):obj.roi(2,1), obj.roi(1,2):obj.roi(2,2), obj.roi(1,3):obj.roi(2,3), obj.roi(1,4):obj.roi(2,4));
            E = gather(E); % converting gpuArray back to normal array
        end
        
        %% Continue to the next iteration. Returns false to indicate that the simulation has terminated
        function state = next(obj, state)
            %% store energy (relative to total energy in source)
            state.diff_energy(state.it) = state.last_step_energy / state.source_energy;
            
            %% check if simulation should terminate
            if (state.diff_energy(state.it) < obj.energy_threshold)
                state.has_next = false;
                state.converged = true;
            elseif (state.it >= state.max_iterations)
                state.has_next = false;
                state.converged = false;
            else
                state.has_next = true;
            end
            
            %% do we need to calculate the last added energy? (only do this once in a while, because of the overhead)
            state.calculate_energy = mod(state.it, obj.energy_calculation_interval)==0;
            
            %% call callback function if neened
            if (mod(state.it, obj.callback_interval)==0 || ~state.has_next) %now and then, call the callback function to give user feedback
                obj.callback(obj, state);
            end
            state.it = state.it+1;
        end
    end
    
    methods(Access = protected)
        %
        % Some helper functions that are not to be called directly
        %
        function [source, source_pos, source_energy] = prepare_source(obj, source, source_pos)
            % prepares the source, returns a cell array with source matrices
            % and a cell array with source position vectors
            % Also verifies that the input is well formed and that
            % all sources fit inside the roi
            %
            if ~iscell(source)
                source = {source};
                source_pos = {source_pos};
            end

            if ~iscell(source_pos) || ~isequal(size(source_pos), size(source))
                error('The number of elements in source and source_pos must match');
            end
            
            % process each source / source_pos combination
            source_energy = 0;
            for c=1:numel(source)             
                % shift source positions so that they are relative to the start of the roi
                pos = simulation.make4(source_pos{c}); %make sure all pos vectors have 4 elements
                source_pos{c} = pos + obj.roi(1,:) - 1;
                if any((source_pos{c} + simulation.make4(size(source{c})) - 1) >  obj.roi(2,:))
                    error('Source does not fit inside the simulation');
                end
                
                % todo: only use non-zero part of source
                % also clear source from GPU memory after usage
                % 
                source{c} = data_array(obj, source{c}); 
                source_energy = source_energy + simulation.energy(source{c});
            end
        end
        
        function d = data_array(obj, data)
            % Creates an array of dimension obj.N. If gpuEnabled is true, the array is created on the gpu
            % Check whether single precision and gpu computation options are enabled
            if obj.single_precision
                p = 'single';
            else
                p = 'double';
            end
            
            if nargin < 2 %no data is specified, fill with zeros
                if obj.gpu_enabled
                    d = zeros(obj.N, p, 'gpuArray');  
                else
                    d = zeros(obj.N, p);
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
            end
        end
        function E = add_sources(state, E, roi, A)
            % roi and source_pos all contain 4-D coordinates (1 for optional polarization)
            % A is an amplitude prefactor to multiply the source by
            
            % process all sources
            for c=1:numel(state.source)
                % determine size of overlap
                pos = state.source_pos{c};
                sz  = simulation.make4(size(state.source{c}));

                %calculate intersection of source and roi
                tlt = max(roi(1,:), pos); %top left corner target
                brt = min(roi(2,:), pos + sz - 1); %bottom right corner target
                if any(tlt > brt)
                    continue; %overlap is empty
                end
                tls = tlt - pos + 1; %top left corner source
                brs = brt - pos + 1; %bottom right corner source
                
                %TODO: the code below is really slow on a gpu because
                % of the indexing. Fortunately we don't call it often in 
                % wavesim, but for PSTD it is _the_ bottleneck
                %
                % it appears that it is already 50% faster if we convert the
                % indices to int64 first! (for a 1x400x1 array set)
                %%
                tls = int64(tls);
                brs = int64(brs);
                tlt = int64(tlt);
                brt = int64(brt);
        
                %add source, multiplied by prefactor A
                E(tlt(1):brt(1), tlt(2):brt(2), tlt(3):brt(3), tlt(4):brt(4)) =...
                        E(tlt(1):brt(1), tlt(2):brt(2), tlt(3):brt(3), tlt(4):brt(4)) + ...
                        A * state.source{c}(tls(1):brs(1), tls(2):brs(2), tls(3):brs(3), tls(4):brs(4));
             end
        end
        
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
            figure(1);
            E = state.E;
            zpos = ceil(size(E, 3)/2);
            imagesc(abs(E(:,:,zpos)));
            title(['Differential energy ' num2str(state.diff_energy(state.it))]);
            drawnow;
        end
        
        %default callback function. Shows real value of field, and total energy evolution
        function default_callback(obj, state)
            figure(1);
            energy = state.diff_energy(1:state.it);
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
            xlabel('y (\lambda / 4)'); ylabel('real(E_x)');
            
            %disp(['Added energy ', num2str(energy(end))]);
            drawnow;
        end
    end
end
classdef simulation
    %Base class for a 2-D wave simulation
    % Ivo M. Vellekoop 2015
    
    properties
        grid; %simgrid object
        roi; %position of simulation area with respect to padded array
        x_range; %x-coordinates of roi
        y_range; %y-coordinates of roi
        z_range; %z-coordinates of roi
        lambda = 1; %wavelength (default = 1 unit)      
        differential_mode = false; %when set to 'true', only the differential field for each iteration is calculated: the fields are not added to get a solution to the wave equation (used for debugging)
        gpu_enabled = false; % flag to determine if simulation are run on the GPU (default: false)
        singlePrecision = false; % flag to determine if single precision is used (default: double)
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

            %copy values from 'options'to 'obj' (only if the field is present in obj)
            obj_fields = fields(obj);
            opt_fields = fields(options);
            fs = intersect(obj_fields, opt_fields);
            for f=1:length(fs)
                fieldname = fs{f};
                obj.(fieldname) = options.(fieldname);
            end
            
            obj.grid = sample.grid;
            obj.roi  = sample.roi;
            obj.x_range = sample.grid.x_range(obj.roi{2});
            obj.x_range = obj.x_range - obj.x_range(1);
            obj.y_range = sample.grid.y_range(obj.roi{1});
            obj.y_range = obj.y_range - obj.y_range(1);
            obj.z_range = sample.grid.z_range(obj.roi{3});
            obj.z_range = obj.z_range - obj.z_range(1);
            if (obj.max_cycles == 0) %default: 1.5 pass
                LN = sqrt(length(obj.x_range).^2 + length(obj.y_range)^2 + length(obj.z_range)^2);
                obj.max_cycles = LN * obj.grid.dx / obj.lambda * 2;
            end
        end
        
        function [E, state] = exec(obj, source, source_pos)
            %% Executes the simulation with the given source
            % The size of 'source' should not exceed the size of the
            % refractive index map, but it can be smaller. In that case
            % the coordinate vector source_pos can be used to indicate
            % the indices of the source within the refractive index map
            % by default, 'source_pos' = [1,1,1]
            tic;
            
            if nargin < 3
                source_pos = [1,1,1];
            end
            source_size = size(source);
            if numel(source_size) < 3
                source_size = [1, source_size];
            end
            
            source_pos = source_pos + [obj.roi{1}(1), obj.roi{2}(1), obj.roi{3}(1)] - [1,1,1];
            source_range = cell(3,1);
            source_range{1} = source_pos(1) + (0:source_size(1)-1);
            source_range{2} = source_pos(2) + (0:source_size(2)-1);
            source_range{3} = source_pos(3) + (0:source_size(3)-1);
            
            if any((source_pos + source_size - [1,1,1]) >  [obj.roi{1}(end), obj.roi{2}(end), obj.roi{3}(end)])
                error('Source does not fit inside the simulation');
            end
            
            %%% prepare state (contains all data that is unique to a single 
            % run of the simulation)
            %
            state.source = data_array(obj, source); %note that source will be converted to 2d automatically if the last dimension is a singleton
            state.source_range = source_range;
            state.it = 1; %iteration
            state.max_iterations = ceil(obj.max_cycles * obj.iterations_per_cycle);
            disp(state.max_iterations);
            state.diff_energy = zeros(state.max_iterations,1);
            state.threshold = obj.energy_threshold * simulation.energy(source); %energy_threshold = fraction of total input energy
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
            end
            E = state.E(obj.roi{1}, obj.roi{2}, obj.roi{3}); %% return only part inside roi. Array remains on the gpu if gpuEnabled = true
        end;
        
        %% Creates an array of dimension obj.grid.N. If gpuEnabled is true, the array is created on the gpu
        function d = data_array(obj, data)
            %% Check whether single precision and gpu computation options are enabled
            if nargin < 2
                if obj.singlePrecision
                    d = zeros(obj.grid.N,'single');
                else
                    d = zeros(obj.grid.N,'double');
                end
            else
                if obj.singlePrecision
                    d = full(single(data));
                else
                    d = full(double(data));
                end
            end
            if obj.gpu_enabled
                d = gpuArray(d);
            end
            if ismatrix(d)
                d = reshape(d, [1, size(d)]); %make 3-D. Note: the result will still be 2-D when the last dimension is a singleton!
            end
        end;
        
        %% Continue to the next iteration. Returns false to indicate that the simulation has terminated
        function state = next(obj, state)
            %% store energy
            state.diff_energy(state.it) = gather(state.last_step_energy);           
            
            %% check if simulation should terminate
            if (state.last_step_energy < state.threshold)
                state.has_next = false;
                state.converged = true;
            elseif (state.it >= state.max_iterations)
                state.has_next = false;
                state.converged = false;
            else
                state.has_next = true;
            end;
            
            %% do we need to calculate the last added energy? (only do this once in a while, because of the overhead)
            state.calculate_energy = mod(state.it, obj.energy_calculation_interval)==0;
            
            %% call callback function if neened
            if (mod(state.it, obj.callback_interval)==0 || ~state.has_next) %now and then, call the callback function to give user feedback
                obj.callback(obj, state);
            end
            state.it = state.it+1;
            
            %% For the differential mode, shut down the source after iteration 1
            if (obj.differential_mode)
                state.source = 0;
            end
        end
    end
    
    methods(Static)
        function en = energy(E_x)
            en= sum(abs(E_x(:)).^2);
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
            threshold = state.threshold;
            E = state.E;
            subplot(2,1,1); plot(1:length(energy),log10(energy),'b',[1,length(energy)],log10(threshold)*ones(1,2),'--r');
            title(length(energy));  xlabel('# iterations'); ylabel('log_{10}(energy added)');
            
            % plot midline along longest dimension
            xpos = ceil(size(E, 2)/2);
            ypos = ceil(size(E, 1)/2);
            zpos = ceil(size(E, 3)/2);
            if xpos > max(ypos, zpos)
                sig = log(abs(E(ypos, :, zpos)));
                sroi = obj.roi{2};
            elseif ypos > max(xpos, zpos)
                sig = log(abs(E(:, xpos, zpos)));
                sroi = obj.roi{1};
            else
                sig = log(abs(E(ypos, xpos, :)));
                sroi = obj.roi{3};
            end
            subplot(2,1,2); plot(1:length(sig), squeeze(sig), sroi(1)*ones(1,2), [min(sig),max(sig)], sroi(end)*ones(1,2), [min(sig),max(sig)]);
            title('midline cross-section')
            xlabel('y (\lambda / 4)'); ylabel('real(E_x)');
            
            %disp(['Added energy ', num2str(energy(end))]);
            drawnow;
        end
        
        function E = add_source(state, E)
            % helper function to add the source to 'E' at the proper
            % offset.
            E(state.source_range{1}, state.source_range{2}, state.source_range{3}) = E(state.source_range{1}, state.source_range{2}, state.source_range{3}) + state.source;
        end
    end
end
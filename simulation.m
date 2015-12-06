classdef simulation
    %Base class for a 2-D wave simulation
    % Ivo M. Vellekoop 2015
    
    properties
        grid; %simgrid object
        roi; %position of simulation area with respect to padded array
        x_range; %x-coordinates of roi
        y_range; %y-coordinates of roi
        lambda = 1; %wavelength (default = 1 unit)        
        gpuEnabled = false; % flag to determine if simulation are run on the GPU (default: false)
        callback = @simulation.default_callback; %callback function that is called for showing the progress of the simulation. Default shows image of the absolute value of the field.
        callback_interval = 50; %the callback is called every 'callback_interval' steps. Default = 5
        max_iterations = 1E4; %1E4; %or when 'max_iterations' is reached. Default 10000
    end
    
    methods
        function obj=simulation(sample, options)
            %% Constructs a simulation object
            %	sample = SampleMedium object
            %   options = simulation options. 

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
        end
        
        function [E, state] = exec(obj, source)
            tic;
            
            %%% prepare state (contains all data that is unique to a single 
            % run of the simulation)
            %
            state.it = 1; %iteration
            state.diff_energy = zeros(obj.max_iterations,1);
            state.threshold = obj.energy_threshold * simulation.energy(source); %energy_threshold = fraction of total input energy
            state.last_step_energy = inf;
            state.calculate_energy = true;
            state.has_next = true;
            
            %% increase source array (which currently the size of the roi)
            % to the full grid size (including boundary conditions)
            if (issparse(source) && ~obj.gpuEnabled)
                state.source = sparse(obj.grid.N(1), obj.grid.N(2));
            else
                state.source = data_array();
            end;
            state.source(obj.roi{1}, obj.roi{2}) = source;

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
            
            E = state.E(obj.roi{1}, obj.roi{2}); %% return only part inside roi. Array remains on the gpu if gpuEnabled = true
        end;
        
        %% Creates an array of dimension obj.grid.N. If gpuEnabled is true, the array is created on the cpu
        function d = data_array(obj)
            %% Check whether gpu computation option is enabled
            if obj.gpuEnabled
                d = gpuArray(obj.grid.N);
            else
                d = zeros(obj.grid.N);
            end
        end;
        
        %% Continue to the next iteration. Returns false to indicate that the simulation has terminated
        function state = next(obj, state)
            %% store energy
            state.diff_energy(state.it) = state.last_step_energy;
            
            %% check if simulation should terminate
            state.it = state.it+1;
            if (state.last_step_energy < state.threshold)
                state.has_next = false;
                state.converged = true;
            elseif (state.it > obj.max_iterations)
                state.has_next = false;
                state.converged = false;
            else
                state.has_next = true;
            end;
            
            %% call callback function if neened
            if (mod(state.it, obj.callback_interval)==0) %now and then, call the callback function to give user feedback
                obj.callback(obj, state);
            end
            
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
        
        %default callback function. Shows real value of field, and total energy evolution
        function default_callback(obj, state)
            figure(1);
            energy = state.diff_energy(1:state.it);
            threshold = state.threshold;
            E = state.E;
            subplot(2,1,1); plot(1:length(energy),log10(energy),'b',[1,length(energy)],log10(threshold)*ones(1,2),'--r');
            title(length(energy));  xlabel('# iterations'); ylabel('log_{10}(energy added)');
            
            sig = log(abs(E(end/2,:)));
            subplot(2,1,2); plot(1:length(sig), sig, obj.roi{2}(1)*ones(1,2), [min(sig),max(sig)], obj.roi{2}(end)*ones(1,2), [min(sig),max(sig)]);
            title('midline cross-section')
            xlabel('y (\lambda / 4)'); ylabel('real(E_x)');
            
            %disp(['Added energy ', num2str(energy(end))]);
            drawnow;
        end
    end
end
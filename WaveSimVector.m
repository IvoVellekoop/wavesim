classdef WaveSimVector < WaveSimBase
    % Implementation of the modified Born series approach for
    % vector waves
    %
    % The only difference that this function implements is in the
    % propagation step
    %
    % Ivo Vellekoop & Gerwin Osnabrugge 2016-2020
    properties
    end
    methods
        function obj=WaveSimVector(refractive_index, options)
            obj@WaveSimBase(refractive_index, options);
            if ~isfield(options,'roi')
                obj.roi(2,4) = 3;  %by default return all 3 polarization planes
            end
            obj.N(4) = 3;
        end
        function E = propagate(obj, E, wiggle)
            % calculate (I - p p^T / (k_0^2+i epsilon)
            % TODO: replace ifft by fft somehow, because ifftn is slow!?
            % apply wiggle phase ramps
            if obj.gpu_enabled
                E = arrayfun(@f_wiggle, E, wiggle.gx, wiggle.gy, wiggle.gz);
            else
                E = f_wiggle(E, wiggle.gx, wiggle.gy, wiggle.gz);
            end
            
            % Fourier transform all vector components separately
            fEx = fftn(E(:,:,:,1));
            fEy = fftn(E(:,:,:,2));
            fEz = fftn(E(:,:,:,3));
            clear E;
            
            % apply propagation kernel
            if obj.gpu_enabled
                [fEx, fEy, fEz] = arrayfun(@f_g0_vector, fEx, fEy, fEz, wiggle.pxe, wiggle.pye, wiggle.pze, obj.k02e);
            else
                [fEx, fEy, fEz] = f_g0_vector(fEx, fEy, fEz, wiggle.pxe, wiggle.pye, wiggle.pze, obj.k02e);
            end
            
            % inverse Fourier transform all vector components
            E = cat(4, ifftn(fEx), ifftn(fEy), ifftn(fEz));
            
            % reverse wiggle phase ramp
            if obj.gpu_enabled
                E = arrayfun(@f_wiggle, E, conj(wiggle.gx), conj(wiggle.gy), conj(wiggle.gz));
            else
                E = f_wiggle(E, conj(wiggle.gx), conj(wiggle.gy), conj(wiggle.gz));
            end
        end
        
        function state = run_algorithm(obj, state)
            
            if(obj.gpu_enabled && obj.usemex && obj.N(4) == 3 )
                
                if(state.max_iterations > 0 && state.max_iterations < inf)
                   maxIter = state.max_iterations;
                else
                   maxIter = 1000;                 
                end
                
                callbackString = 'NoCallback';
                if(strcmp(func2str(obj.callback), 'Simulation.energy_added_disp_callback'))
                    callbackString = 'EnergyAddedDisp';
                elseif(strcmp(func2str(obj.callback), 'Simulation.energy_added_callback'))
                    callbackString = 'EnergyAdded';             
                end
                
                k02eSingle = single(gather(obj.k02e));
                epsilonSingle = single(gather(obj.epsilon));
                Nwiggles = numel(obj.wiggles);
                maxIter = double(gather(maxIter));

                [state.E, state.diff_energy] = RunWaveSimAlgorithmMex(obj.gamma, Nwiggles, obj.wiggles, state.source, ...
                epsilonSingle, k02eSingle, obj.roi, maxIter, obj.energy_threshold, callbackString, obj.callback_interval);
                
                %% Final computional steps
                % divide field by gamma to convert E' -> E
                state.E = state.E ./ obj.gamma; 

                % crop field to match roi size
                state.E = obj.crop_field(state.E);
                state.it = length(state.diff_energy) +1;
                state.converged = state.it < maxIter;
            else 
                state = run_algorithm@WaveSimBase(obj, state);
            end
        end
    end
end

classdef WaveSim < WaveSimBase
    % Implementation of the modified Born series approach for
    % scalar waves
    %
    % The only difference that this function implements is in the
    % propagation step
    %
    % Ivo Vellekoop & Gerwin Osnabrugge 2016-2020
    properties
    end
    methods
        function obj = WaveSim(refractive_index, options)
            obj@WaveSimBase(refractive_index, options);
        end
        function E = propagate(obj, E, wiggle)
            % calculate (I - p p^T / (k_0^2+i epsilon)
            
            % prepare coordinate vectors (shifted by quarter of a pixel if
            % wiggling)
            px2 = wiggle.pxe.^2;
            py2 = wiggle.pye.^2;
            pz2k = wiggle.pze.^2 - obj.k02e;
                        
            % perform modified Fourier transform (applies wiggle phase ramps)
            E = obj.mfft(E, wiggle);
            
            % apply propagation kernel
            if obj.gpu_enabled
                E = arrayfun(@f_g0_scalar, E, px2, py2, pz2k);
            else
                E = f_g0_scalar(E, px2, py2, pz2k);
            end
            
            % Perform modified inverse Fourier transform (reverses wiggle phase ramps)
            E = obj.mifft(E, wiggle);    
        end
        % modified Fourier transforms (used for acyclic convolution)
        function E = mfft(obj, E, wig)
            % Modified Fourier transform: additionally applies wiggle phase
            % ramps to field. Maps field from real space to a shifted
            % Fourier space 
            if obj.gpu_enabled
                E = arrayfun(@f_wiggle, E, wig.gx, wig.gy, wig.gz);
            else
                E = f_wiggle(E, wig.gx, wig.gy, wig.gz);
            end
            E = fftn(E);
        end       
        function E = mifft(obj, E, wig)
            % Modified inverse Fourier transform: additionally reverses wiggle
            % phase ramps applied to field. Maps field from a shifted Fourier 
            % space to real space  
            E = ifftn(E);
            if obj.gpu_enabled
                E = arrayfun(@f_wiggle, E, conj(wig.gx), conj(wig.gy), conj(wig.gz));
            else
                E = f_wiggle(E, conj(wig.gx), conj(wig.gy), conj(wig.gz));
            end
        end

        function state = run_algorithm(obj, state)
            
            if(obj.gpu_enabled && obj.usemex )
                if(state.max_iterations > 0 && state.max_iterations < inf)
                   maxIter = state.max_iterations;
                else
                   maxIter = 100000;                 
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

                [state.E, state.diff_energy] = RunWaveSimAlgorithmMex(false, obj.gamma, Nwiggles, obj.wiggles, state.source, ...
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


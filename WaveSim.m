classdef WaveSim < WaveSimBase
    % Implementation of the modified Born series approach for
    % scalar waves
    %
    % The only difference that this function implements is in the
    % propagation step
    %
    % Ivo M. Vellekoop 2018
    properties
    end
    methods
        function obj = WaveSim(sample, options)
            obj@WaveSimBase(sample, options);
        end
        function E = propagate(obj, E)
            % calculate (I - p p^T / (k_0^2+i epsilon)
            if obj.gpu_enabled
                E = ifftn(arrayfun(f_g0_scalar, fftn(E), obj.pxe.^2, obj.pye.^2, obj.pze.^2 - obj.k02e));
            else
                E = ifftn(f_g0_scalar(fftn(E), obj.pxe.^2, obj.pye.^2, obj.pze.^2 - obj.k02e));
            end    
        end
    end
end


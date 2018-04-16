classdef WaveSimVector < WaveSimBase
    % Implementation of the modified Born series approach for
    % vector waves
    %
    % The only difference that this function implements is in the
    % propagation step
    %
    % Ivo M. Vellekoop 2018
    properties
    end
    methods
        function obj=WaveSimVector(sample, options)
            obj@WaveSimBase(sample, options);
            obj.roi(2,4) = 3; %3-polarization planes
            obj.N(4)     = 3;
        end
        function E = propagate(obj, E)
            % calculate (I - p p^T / (k_0^2+i epsilon)
            fEx = fftn(E(:,:,:,1));
            fEy = fftn(E(:,:,:,2));
            fEz = fftn(E(:,:,:,3));
            if obj.gpu_enabled
                [fEx, fEy, fEz] = arrayfun(f_g0_vector, fEx, fEy, fEz, obj.pxe, obj.pye, obj.pze, obj.k02e);
            else
                [fEx, fEy, fEz] = f_g0_vector(fEx, fEy, fEz, obj.pxe, obj.pye, obj.pze, obj.k02e);
            end    
            E(:,:,:,1) = ifftn(fEx);
            E(:,:,:,2) = ifftn(fEy);
            E(:,:,:,3) = ifftn(fEz);
        end
    end
end

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
            % TODO: replace fEx with E directly to save memory
            % apply wiggle phase ramps

            % Apply (modified) Fourier transform to all vector components 
            % separately (also applies wiggle phase ramps)
            fEx = obj.mfft(E(:,:,:,1),wiggle);
            fEy = obj.mfft(E(:,:,:,2),wiggle);
            fEz = obj.mfft(E(:,:,:,3),wiggle);
            
            % apply propagation kernel
            if obj.gpu_enabled
                [fEx, fEy, fEz] = arrayfun(@f_g0_vector, fEx, fEy, fEz, wiggle.pxe, wiggle.pye, wiggle.pze, obj.k02e);
            else
                [fEx, fEy, fEz] = f_g0_vector(fEx, fEy, fEz, wiggle.pxe, wiggle.pye, wiggle.pze, obj.k02e);
            end
            
            % apply (modified) inverse Fourier transform to all vector components
            E(:,:,:,1) = obj.mifft(fEx,wiggle);
            E(:,:,:,2) = obj.mifft(fEy,wiggle);
            E(:,:,:,3) = obj.mifft(fEz,wiggle);
        end
    end
end

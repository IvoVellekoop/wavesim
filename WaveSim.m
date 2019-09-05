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
        function E = propagate(obj, E, wiggle)
            % calculate (I - p p^T / (k_0^2+i epsilon)
            
            % prepare coordinate vectors (shifted by quarter of a pixel if
            % wiggling)
            px2 = wiggle.pxe.^2;
            py2 = wiggle.pye.^2;
            pz2k = wiggle.pze.^2 - obj.k02e;
            
            % perform modified Fourier transform (applies wiggle phase ramps)
            E = obj.fft_wiggle(E, wiggle);
            
            % apply propagation kernel
            if obj.gpu_enabled
                E = arrayfun(@f_g0_scalar, E, px2, py2, pz2k);
            else
                E = f_g0_scalar(E, px2, py2, pz2k);
            end
            
            % Perform modified inverse Fourier transform (reverses wiggle phase ramps)
            E = obj.ifft_wiggle(E, wiggle);          
        end
    end
end


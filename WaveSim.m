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
            E = fft_wiggle(E, wiggle, obj.gpu_enabled);
            
            % apply propagation kernel
            if obj.gpu_enabled
                E = arrayfun(@f_g0_scalar, E, px2, py2, pz2k);
            else
                E = f_g0_scalar(E, px2, py2, pz2k);
            end
            
            % Perform modified inverse Fourier transform (reverses wiggle phase ramps)
            E = ifft_wiggle(E, wiggle,obj.gpu_enabled);          
%             if obj.gpu_enabled
%                 % wiggle, then Fourier transform
%                 E = fftn(arrayfun(@f_wiggle, E, wiggle.gx, wiggle.gy, wiggle.gz));
%                 
% 
%                 % apply propagation kernel
%                 E = arrayfun(@f_g0_scalar, E, px2, py2, pz2k);
%                 
%                 % Fourier tranform back, then wiggle back
%                 E = arrayfun(@f_wiggle, ifftn(E), conj(wiggle.gx), conj(wiggle.gy), conj(wiggle.gz));
%             else
%                 E = fftn(f_wiggle(E, wiggle.gx, wiggle.gy, wiggle.gz));
%                 E = f_g0_scalar(E, px2, py2, pz2k);
%                 E = f_wiggle(ifftn(E), conj(wiggle.gx), conj(wiggle.gy), conj(wiggle.gz));
%             end
        end
    end
end


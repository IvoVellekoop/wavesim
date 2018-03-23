classdef wavesimv < wavesim_base
    %Simulation of the 2-D or 3-D Maxwell's equations wave equation using a Born series approach
    %Only valid for non-magnetic media.
    % Ivo M. Vellekoop 2017
    properties
    end
    methods
        function obj=wavesimv(sample, options)
            obj@wavesim_base(sample, options);
            obj.roi(4,2) = 3; %3-polarization planes
            obj.N(4)     = 3;
            % Set up propagation function (convolution with Green's function)
            % the Green's function is pre-multiplied by epsilon to reduce the
            % number of computations a bit
            %
            k0e = obj.data_array(obj.k^2 + 1.0i * obj.epsilon);
            px = obj.data_array(obj.grid.px_range);
            py = obj.data_array(obj.grid.py_range);
            pz = obj.data_array(obj.grid.pz_range);
            eps = obj.data_array(obj.epsilon);
            
            if obj.gpu_enabled
                obj.propagate = @(E) wavesimv.propagate_gpu(E, px, py, pz, k0e, eps);
            else
%                p = zeros(obj.N);
%                [p(:,:,:,1), p(:,:,:,2), p(:,:,:,3)] = ndgrid(obj.px, obj.py, obj.pz);
                g0_k = obj.epsilon ./ (px.^2+py.^2+pz.^2-k0e);
                obj.propagate = @(E) wavesimv.propagate_cpu(E, px, py, pz, k0e, g0_k);
            end
        end
    end
    methods (Static)
        function E = propagate_cpu(E, px, py, pz, k0e, g0_k)
            % calculate (I - p p^T / (k_0^2+i epsilon)
            fEx = fftn(E(:,:,:,1));
            fEy = fftn(E(:,:,:,2));
            fEz = fftn(E(:,:,:,3));
            div = (px .* fEx + py .* fEy + pz .* fEz) / k0e; %divergence term
            fEx = g0_k .* (fEx - px .* div);
            fEy = g0_k .* (fEy - py .* div);
            fEz = g0_k .* (fEz - pz .* div);
            E(:,:,:,1) = ifftn(fEx);
            E(:,:,:,2) = ifftn(fEy);
            E(:,:,:,3) = ifftn(fEz);
        end
        function E = propagate_gpu(E, px, py, pz, k0e, epsilon)
            % calculate (I - p p^T / (k_0^2+i epsilon)
            fEx = fftn(E(:,:,:,1));
            fEy = fftn(E(:,:,:,2));
            fEz = fftn(E(:,:,:,3));
            [fEx, fEy, fEz] = arrayfun(@f_g0, fEx, fEy, fEz, px, py, pz, k0e, epsilon);
            E(:,:,:,1) = ifftn(fEx);
            E(:,:,:,2) = ifftn(fEy);
            E(:,:,:,3) = ifftn(fEz);
        end
    end
end

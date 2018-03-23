classdef wavesimv < wavesim_base
    %Simulation of the 2-D or 3-D Maxwell's equations wave equation using a Born series approach
    %Only valid for non-magnetic media.
    % Ivo M. Vellekoop 2017
    properties
        px;
        py;
        pz;
        g0_k;
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
            k0e = obj.k^2 + 1.0i * obj.epsilon;
            obj.px = obj.grid.px_range/sqrt(k0e);
            obj.py = obj.grid.py_range/sqrt(k0e);
            obj.pz = obj.grid.pz_range/sqrt(k0e);
            obj.g0_k = 1./(k0e/obj.epsilon * (obj.px.^2+obj.py.^2+obj.pz.^2-1));
            obj.propagate = @(E) wavesimv.propagate(obj, E);
        end
    end
    methods (Static)
        function E = propagate(obj, E)
            % calculate (I - p p^T / (k_0^2+i epsilon)
            fEx = fftn(E(:,:,:,1));
            fEy = fftn(E(:,:,:,2));
            fEz = fftn(E(:,:,:,3));
            div = obj.px .* fEx + obj.py .* fEy + obj.pz .* fEz; %divergence term
            fEx = obj.g0_k .* (fEx - obj.px .* div);
            fEy = obj.g0_k .* (fEy - obj.py .* div);
            fEz = obj.g0_k .* (fEz - obj.pz .* div);
            E(:,:,:,1) = ifftn(fEx);
            E(:,:,:,2) = ifftn(fEy);
            E(:,:,:,3) = ifftn(fEz);
        end
    end
end
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
            
            % Set up propagation function (convolution with Green's function)
            % the Green's function is pre-multiplied by epsilon to reduce the
            % number of computations a bit
            %
            px2 = obj.grid.px_range.^2/obj.epsilon;
            py2 = obj.grid.py_range.^2/obj.epsilon;
            pz2e = obj.grid.pz_range.^2/obj.epsilon - (obj.k^2 + 1.0i * obj.epsilon)/obj.epsilon;
            if obj.gpu_enabled
                % generate g0_k on the GPU (saves memory and memory
                % bandwidth, making it faster than pre-calculating g0_k)
                multiply_g0_k = @(E, px2, py2, pz2e) ...
                    E ./(px2+py2+pz2e);
                obj.propagate = @(E) ...
                    ifftn(arrayfun(multiply_g0_k, fftn(E), px2, py2, pz2e));
            else
                % pre-calculate g0_k when using the CPU (unfortunately this 
                % distinction is needed because arrayfun for the CPU does not 
                % support singleton expansion like arrayfun for the GPU!)
                % also, arrayfun is incredibly slow on the CPU!
                g0_k = 1./(px2+py2+pz2e);
                obj.propagate = @(E) ...
                    ifftn(g0_k .* fftn(E));
            end
        end
    end
end


function E = f_g0_scalar(E, pxe2, pye2, pze2k0e)
% Propagation function (convolution with Green's function)
% the Green's function is pre-multiplied by epsilon to reduce the
% number of computations a bit
%
    E = E ./(pxe2+pye2+pze2k0e);
end

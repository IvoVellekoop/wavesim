%RunWaveSimAlgorithmMex
%   [E, diff_energy] = RunWaveSimAlgorithmMex(gamma, wiggles, source, ...
%                epsilon, k02e, roi, max_iter, energy_threshold);
%   Calculates the energy using the modified Born series approach for
%   vector waves, using the wiggle phase ramps (gx, gy, gz) and shifted 
%   Fourier corordinates (pxe, pye, pze) for the anti-cyclic convolution. 
%   Energy is added iteratively.
%
%   Output:
%   E [Ny x Nx x Nz x 3] - complex gpu array
%   diff_energy - the total added energy per iteration
%
%   Input:
%   gamma [Ny x Nx x Nz] - complex gpu array calculated from Medium
%   refractive index map
%   wiggles [8 x 1] - cell array containing wiggle structs with elements
%       gx [1 x Nx] complex gpu array
%       gy [Ny x 1] complex gpu array
%       gx [1 x 1 x Nz] complex gpu array
%       pxe [1 x Nx] real gpu array
%       pye [Ny x 1] real gpu array
%       pze [1 x 1 x Nz] real gpu array
%   source - Source class
%   epsilon - real scalar in single precision
%   k02e - complex scalar in single precision
%   roi - [2x4] double matrix containing start and end index of the
%   simulation roi
%   max_iter - (double) maximum number of iterations (over a whole wiggle sequence of length 8)
%   energy_threshold - (double) the threshold for the factor diff_energy(end)/diff_energy(1)
%   signaling convergence

%   MEX-File function.
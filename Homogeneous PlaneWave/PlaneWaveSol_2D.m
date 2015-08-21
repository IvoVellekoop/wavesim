function [ sol ] = PlaneWaveSol_2D( k, theta, x, y )
% ANALYTIC_SOL_2D 
% Compute the analytic solution to the Helmholtz equation in 2D in
% homogeneous medim
% param:
%   k: wavenumber k=k0*n (d²u/dx² + k² * u = 0)
%   Theta: parameter of the solution.
%   x, y: set of the points

    k1 = k * cos(theta);
    k2 = k * sin(theta);
    sol = exp(1.0i * (k1 * x + k2 * y + pi/2));
    
end


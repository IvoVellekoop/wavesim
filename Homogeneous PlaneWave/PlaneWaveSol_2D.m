function [ sol ] = PlaneWaveSol_2D( k, theta, x, y )
% ANALYTIC_SOL_2D 
% Compute the analytic solution to the Helmholtz equation in 2D in
% homogeneous medim
% param:
%   k: wavenumber k=k0*n (d²u/dx² + k² * u = 0)
%   Theta: parameter of the solution.
%   x, y: set of the points

    kx = k * cos(theta);
    ky = k * sin(theta);
    sol = exp(1.0i * (kx * x + ky * y -pi/2 ));%
    
end


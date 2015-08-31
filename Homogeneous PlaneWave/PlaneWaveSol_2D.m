function [ sol ] = PlaneWaveSol_2D( k, theta,initalphase,h, x, y )
% ANALYTIC_SOL_2D 
% Compute the analytic solution to the Helmholtz equation in 2D in
% homogeneous medim
% param:
%   k: wavenumber k=k0*n (d²u/dx² + k² * u = 0)
%   Theta: parameter of the solution.
%   x, y: set of the points

    kx = k * cos(theta);
    ky = k * sin(theta);
    %sol = exp(1.0i * (kx * x + ky * y + pi/2));
    %sol = ( - exp(1.0i * abs(kx * x + ky * y - initalphase)) + exp(1.0i*abs(kx * x + ky * y))); % normalized
    %sol = (1/(2*k^2))*( - exp(1.0i * abs(kx * x + ky * y - initalphase)) + exp(1.0i*abs(kx * x + ky * y))); % [o h]
    %sol = (1/(2*k^2))*( - exp(1.0i * (kx * x + ky * y - initalphase/2)) + exp(1.0i*(kx * x + ky * y + initalphase/2))); % [-h/2 h/2]
    sol = (-2.0i/(k^3*h))*(exp(1.0i * (kx * x + ky * y - initalphase/2)) + exp(1.0i*(kx * x + ky * y + initalphase/2))- 2*exp(1.0i*(kx * x + ky * y))); % Triangle source [-h/2 h/2]
end


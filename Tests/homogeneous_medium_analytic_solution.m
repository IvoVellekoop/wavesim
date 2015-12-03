function [ sol ] = homogeneous_medium_analytic_solution( k, theta, h, x, y )
% ANALYTIC_SOL_2D
% Compute the analytic solution to the Helmholtz equation in 2D in
% homogeneous medim
% param:
%   k: wavenumber k=k0*n (d²u/dx² + k² * u = 0)
%   theta: angle of the wave
%   h: grid sixe
%   x, y: meshgrid of coordinates (todo: change to bsxfun)
kx = k * cos(theta);
ky = k * sin(theta);

%finding the exact solution is not trivial because
%the source is not an exact delta function. 
%The best correspondence with the discrete representation
%is found by assuming that the source is a sinc function with a width
%matching the grid spacing (i. e. the source is band-limited)

%todo: find out why factor 1/4pi is correct (solution with Mathematica
%gave 1/8pi prefactor for expint term?)
sol = 1.0i*h/(2*k)*exp(1.0i * (kx * x + ky * y))... %<--propagating plane wave.
    -h/(4*pi*k) * (...
    exp(1.0i * (kx * x + ky * y)) .* (  expint(-1.0i * (-h*kx+pi) * x / h) - expint(1.0i * (h*kx+pi) * x/h)) -...
    exp(-1.0i * (kx * x + ky * y)) .* ( -expint(-1.0i * (h*kx-pi)* x / h)  + expint(-1.0i* (h*kx+pi) * x/h)));
end


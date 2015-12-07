function [ sol ] = homogeneous_medium_analytic_solution( k, h, x )
% ANALYTIC_SOL
% Compute the analytic solution to the Helmholtz equation in 2D in
% homogeneous medim
% param:
%   k: wavenumber k=k0*n (d²u/dx² + k² * u = 0)
%   h: grid sixe
%   x: range of coordinates

%finding the exact solution is not trivial because
%the source is not an exact delta function. 
%The best correspondence with the discrete representation
%is found by assuming that the source is a sinc function with a width
%matching the grid spacing (i. e. the source is band-limited)

%todo: find out why factor 1/4pi is correct (solution with Mathematica
%gave 1/8pi prefactor for expint term?)
phi = k * x;
x(abs(x)<1E-100) = 1E-100; %dirty way to avoid nan
sol = 1.0i*h/(2*k)*exp(1.0i * phi)... %<--propagating plane wave.
    -h/(4*pi*k) * (...
    exp(1.0i * phi) .* (  expint(1.0i * (k-pi/h) * x) - expint(1.0i * (k+pi/h) * x)) -...
    exp(-1.0i * phi) .* ( -expint(-1.0i * (k-pi/h)* x)  + expint(-1.0i* (k+pi/h) * x)));
end


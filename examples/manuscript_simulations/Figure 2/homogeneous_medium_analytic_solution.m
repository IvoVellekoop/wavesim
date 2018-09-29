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
%see Mathematica file for derivation

% To determine peak value at x=0, realize that Ei(x) ~ ln(x) for x close to
% 0 and find E(0) = -h/(4*pi*k)*2*log((pi/h+k)/(pi/h-k))  + i*h/(2*k)
% which can be rewritten as E(0) = -h/(2*pi*k)*log((1+h*k/pi)/(1-h*k/pi)) + i*h/(2*k)
% E(0) = -h/(pi*k) * atanh(h*k/pi) + i*h/(2*k)
% E(0) = h/(pi*k) * (i*pi/2 -atanh(h*k/pi))
% E(0) = i*h/(2*k) * (1 + 2i*atanh(h*k/pi)/pi)
% note: h*k/pi = 2*h/lambda  ==> Sampling rate with respect to Nyquist
% note: the Ei part (sol-plane wave component) is real and rapidly
% oscillating
phi = k * x;
x(abs(x)<1E-100) = 1E-100; %dirty way to avoid nan
sol = 1.0i*h/(2*k)*exp(1.0i * phi)... %<--propagating plane wave.
    -h/(4*pi*k) * (...
    exp(1.0i * phi) .* (  expint(1.0i * (k-pi/h) * x) - expint(1.0i * (k+pi/h) * x)) -...
    exp(-1.0i * phi) .* ( -expint(-1.0i * (k-pi/h)* x)  + expint(-1.0i* (k+pi/h) * x)));
end


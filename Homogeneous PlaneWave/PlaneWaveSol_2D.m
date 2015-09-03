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
    %sol = exp(1.0i * (kx * x + ky * y -pi/2));
    %sol = ( - exp(1.0i * abs(kx * x + ky * y - initalphase)) + exp(1.0i*abs(kx * x + ky * y))); % normalized
    %sol = (1/(2*k^2))*( - exp(1.0i * abs(kx * x + ky * y - initalphase)) + exp(1.0i*abs(kx * x + ky * y))); % [o h]
    %sol = (1/(2*k^2))*( - exp(1.0i * (kx * x + ky * y - initalphase/2)) + exp(1.0i*(kx * x + ky * y + initalphase/2))); % [-h/2 h/2]
    
    % delta source
    %sol = (1/(2*k^2))* (exp(1.0i*(kx * x + ky * y))); 
    
    % Triangle source [-h/2 h/2]
    %sol = (-2.0i/(k^3*h))*(exp(1.0i * (kx * x + ky * y - initalphase/2))...
    %                           + exp(1.0i*(kx * x + ky * y + initalphase/2))- 2*exp(1.0i*(kx * x + ky * y))); 
     sol = ( 8.0i*sin(initalphase/4)^2/(k^3*h) )*exp(1.0i * abs( kx * x + ky * y ));
    
    % Triangle source [-h h]
    %sol = ( (-1.0i*(cos(initalphase)-1))/(k^3*h) )*exp(1.0i * abs( kx * x + ky * y ));
    
    %Sinc function
    %sol = (h/(4*pi*k))*exp(-1.0i * abs(kx * x + ky * y)).*(...
    %            (sgn(x-x')-1)*(-expint(-1.0i*((initalphase-pi)*x')./h)+expint((-1.0i*(initalphase+pi)*x')./h))...
    %            +((sgn(x-x')+1)*exp(2.0i * abs(kx * x + ky * y)))...
    %            .*(-expint(-1.0i*((pi-initalphase)*x')./h)+expint((+1.0i*(initalphase+pi)*x')./h)));
    %            
    %Sinc function (having no abs(x-x'))
    %sol = (h^2/(2*pi^2))*exp(1.0i *(kx * x + ky * y))*(...                             
    %            -((1.0i*pi^2/initalphase)-(2*pi*atan(initalphase/pi)/initalphase))...
    %            +((1.0i*pi^2/initalphase)-(2*pi*atanh(initalphase/pi)/initalphase)));
    %
    %Sinc function (Intergarted from -Inf to 0)
    
    %sol = (h/(2*k*pi))*exp(1.0i * (kx * x + ky * y))*...
    %                ( (-expint(-1.0i*((initalphase-pi)*x')./h)) - expint((+1.0i*(initalphase+pi)*x')./h)+...
    %                        log2(1.0i*(-initalphase+pi)/h) - log2(-1.0i*(initalphase+pi)/h) );
end


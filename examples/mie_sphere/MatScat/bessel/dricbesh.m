function [xi] = dricbesh(nu,k,z)
%DRICBESH First derivative of the Riccati-Bessel function of the third kind.
%   XI = DRICBESH(NU,K,Z), for K = 1 or 2, computes the first derivative 
%   of the Riccati-Bessel function of the third kind for each element of 
%   the complex array Z.
%
%   XI = DRICBESH(NU,Z) uses K = 1.%
%
%   See also:
%   DRICBESJ, DRICBESY, RICBESJ, RICBESY, RICBESH, SBESSELH
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)


if nargin == 2
    z = k;
    k = 1;
end %if nargin == 2

xi = 0.5*sqrt(pi/2/z)*besselh(nu+0.5,k,z) ...
    + sqrt(pi*z/2)*dbesselh(nu+0.5,k,z);

end


function [xi] = ricbesh(nu,k,z)
%RICBESH Riccati-Bessel function of the third kind.
%   XI = RICBESH(NU,K,Z), for K = 1 or 2, computes the Riccati-Bessel 
%   function of the third kind for each element of the complex array Z.
%
%   XI = RICBESH(NU,Z) uses K = 1.
%
%   See also:
%   RICBESJ, RICBESY, DRICBESJ, DRICBESY, DRICBESH, SBESSELH
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

if nargin == 2
    z = k;
    k = 1;
end %if nargin == 2

xi = sqrt(pi*z/2)*besselh(nu+0.5,k,z);

end


function [h] = dbesselh(nu,k,z)
%DBESSELH First derivative of the Bessel function of the third kind.
%   H = DBESSELH(NU,K,Z), for K = 1 or 2, computes the derivative of the
%   Hankel function, dH^(K)_NU(Z)/dZ, for each element of the complex 
%   array Z.
%
%   H = DBESSELH(NU,Z) uses K = 1.
%
%   See also:
%   DBESSELJ, DBESSELY
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

if nargin == 2
    z = k;
    k = 1;
end %if nargin == 2
    
h = 0.5*(besselh(nu-1,k,z) - besselh(nu+1,k,z));

end

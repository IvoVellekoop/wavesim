function [h] = sbesselh(nu,k,z)
%SBESSELH Spherical Bessel function of the third kind.
%
%   H = SBESSELH(NU,K,Z), for K = 1 or 2, computes the spherical Bessel
%   function of the third kind, h^(K)_NU(Z), for each element of the 
%   complex array Z.
%
%   H = SBESSELH(NU,Z) uses K = 1.
%
%   See also:
%   SBESSELJ, SBESSELY, RICBESH, DRICBESH
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

if nargin == 2
    z = k;
    k = 1;
end %if nargin == 2
    
h = ricbesh(nu,k,z)./z;

end

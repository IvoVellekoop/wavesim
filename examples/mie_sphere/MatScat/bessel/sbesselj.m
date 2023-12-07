function [j] = sbesselj(nu,z)
%SBESSELH Spherical Bessel function of the first kind.
%
%   J = SBESSELJ(NU,Z), computes the spherical Bessel function of the first
%   kind, j_NU(Z), for each element of the complex array Z.
%
%   See also:
%   SBESSELH, SBESSELY, RICBESJ, DRICBESJ
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

j = ricbesj(nu,z)./z;

end

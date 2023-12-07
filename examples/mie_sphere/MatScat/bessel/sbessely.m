function [y] = sbessely(nu,z)
%SBESSELH Spherical Bessel function of the second kind.
%
%   Y = SBESSELY(NU,Z), computes the spherical Bessel function of the
%   second kind, y_NU(Z), for each element of the complex array Z.
%
%   See also:
%   SBESSELJ, SBESSELH, RICBESY, DRICBESY
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

y = -ricbesy(nu,z)./z;

end

function [S] = dricbesj(nu,z)
%DRICBESJ First derivative of the Riccati-Bessel function of the first kind.
%   S = DRICBESJ(NU,Z) is the first derivative of the Riccati-Bessel 
%   function of the first kind for each element of the complex array Z.
%
%   See also:
%   DRICBESY, DRICBESH, RICBESJ, RICBESY, RICBESH, SBESSELH
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

S = 0.5*sqrt(pi/2/z)*besselj(nu+0.5,z) + sqrt(pi*z/2)*dbesselj(nu+0.5,z);

end


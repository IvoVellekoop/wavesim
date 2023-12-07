function [C] = ricbesy(nu,z)
%RICBESY Riccati-Bessel function of the second kind.
%   C = RICBESY(NU,Z) is the Riccati-Bessel function of the second kind for
%   each element of the complex array Z.
%
%   See also:
%   RICBESJ, RICBESH, DRICBESJ, DRICBESY, DRICBESH, SBESSELH
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

C = -sqrt(pi*z/2)*bessely(nu+0.5,z);

end


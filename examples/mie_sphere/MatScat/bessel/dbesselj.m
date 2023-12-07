function [j] = dbesselj(nu,z)
%DBESSELJ First derivative of the Bessel function of the first kind.
%   J = DBESSELJ(NU,Z) is the first derivative of Bessel function of the 
%   first kind, dJ_NU(Z)/dZ, for each element of the complex array Z.
%
%   See also:
%   DBESSELY, DBESSELH
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

j = 0.5*(besselj(nu-1,z) - besselj(nu+1,z));

end


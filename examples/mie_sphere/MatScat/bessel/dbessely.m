function [y] = dbessely(nu,z)
%DBESSELY First derivative of the Bessel function of the second kind.
%   Y = DBESSELY(NU,Z) is the first derivative of Bessel function of the 
%   second kind, Y_NU(Z)/dZ, for each element of the complex array Z.
%
%   See also:
%   DBESSELJ, DBESSELH
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

y = 0.5*(bessely(nu-1,z) - bessely(nu+1,z));

end


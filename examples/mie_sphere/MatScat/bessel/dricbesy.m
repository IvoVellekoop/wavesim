function [C] = dricbesy(nu,z)
%DRICBESY First derivative of the Ricatti-Bessel function of the second kind.
%   C = DRICBESY(NU,Z) is the first derivative of the Riccati-Bessel 
%   function of the second kind for each element of the complex array Z.
%
%   See also:
%   DRICBESJ, DRICBESH, RICBESJ, RICBESY, RICBESH, SBESSELH
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

C = -0.5*sqrt(pi/2/z)*bessely(nu+0.5,z) - sqrt(pi*z/2)*dbessely(nu+0.5,z);

end


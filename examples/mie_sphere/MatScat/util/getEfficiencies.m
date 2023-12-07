function [Q] = getEfficiencies( C, r, dim )
%GETEFFICIENCIES Calculates the efficiencies for given cross sections
%
%   [Q] = GETEFFICIENCIES(C,R,DIM) calculates the extinction, scattering
%   and absorption efficiencies Q with given cross sections C for a
%   DIM-dimensional sphere with radius R.
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

if dim == 2
    Cgeo = 2*r;
elseif dim == 3
    Cgeo = pi*r^2;
end %if dim

Q.ext = C.ext/Cgeo;
Q.sca = C.sca/Cgeo;
Q.abs = C.abs/Cgeo;
end

function [ theta, phi, r ] = xcart2sph( x, y, z )
%XCART2SPH Transform cartesian to spherical coordinates using a different
%   convention than the standard cart2sph function.
%
%   [THETA,PHI,R] = XCART2SPH(X,Y,Z) transforms corresponding elements of
%   data stored in Cartesian coordinates X,Y,Z to spherical coordinates
%   coordinates (THETA,PHI,R). The arrays X,Y, and Z must be the same size 
%   (or any of them can be scalar). THETA and PHI are returned in radians.
%
%   See also:
%   CART2SPH
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

theta = pi/2. - atan(z./sqrt(x.^2 + y.^2));
phi = atan2(y,x);
r = sqrt(x.^2 + y.^2 + z.^2);

end


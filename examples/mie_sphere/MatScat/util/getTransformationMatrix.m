function [R] = getTransformationMatrix( phi, theta )
%GETTRANSFORMATIONMATRIX Returns the matrix to transform a 3D vector from
%   spherical/cylindrical coordinate basis to cartesian coordinate basis.
%
%   [R] = NFCYL(PHI,THETA) returns the 3x3 transformation matrix to express
%   a 3D vector in cartesian coordinate basis. If PHI and THETA is given
%   the original vector basis is assumed to be in spherical coordinates, if
%   only PHI is submitted original vector basis is in cylindrical
%   coordinates.
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

R = zeros(3,3);

if ~exist('theta', 'var')
    R(1,1) = cos(phi);
    R(1,2) = -sin(phi);
    R(1,3) = 0.;
    R(2,1) = sin(phi);
    R(2,2) = cos(phi);
    R(2,3) = 0.;
    R(3,1) = 0.;
    R(3,2) = 0.;
    R(3,3) = 1.;
else %~exist('theta', 'var')
    R(1,1) = sin(theta)*cos(phi);
    R(1,2) = cos(theta)*cos(phi);
    R(1,3) = -sin(phi);
    R(2,1) = sin(theta)*sin(phi);
    R(2,2) = cos(theta)*sin(phi);
    R(2,3) = cos(phi);
    R(3,1) = cos(theta);
    R(3,2) = -sin(theta);
    R(3,3) = 0;
end %if ~exist('theta', 'var')

end


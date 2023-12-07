function [T, C, ang] = calccyl_multi( r, ns, nm, lambda, xc, yc, nang, ...
    iphi, varargin )
%CALCCYL_MULTI Calculates the amplitude scattering matrix and cross 
%   sections for the scattering of electromagnetic radiation by multiple 
%   infinite circular cylinders at perpendicular incidence.
%
%   [T,C,ANG] = CALCCYL_MULTI(R,NS,NM,LAMBDA,XC,YC,NANG,IPHI,VARARGIN) 
%   calculates the amplitude scattering matrix T and cross sections C for 
%   the scattering of a monochromatic electromagnetic wave with vacuum 
%   wavelength LAMBDA at perpendicular incidence by many cylinders with 
%   radius R and (complex) refractive index NS. NM is the refractive index 
%   of the surrounding medium. The coordinates of the cylinders are given 
%   in XC and YC, the cylinders may not overlap. The coordinate lists XC 
%   and YC have to be the same size. The incident light is propagating in 
%   the xy-plane at an angle IPHI (in degree) according to the positive 
%   x-axis. The cylinder axes are aligned along the z-axis.
%
%   VARARGIN can take the following keywords:
%    'HankelFunction'       : HK - Use hankel function H_n^(HK) in
%       calculations (default HK = 1).
%       WARNING: The complex refractive index is defined as:
%       n = nr + ni*i (HK == 1)
%       n = nr - ni*i (HK == 2)
%    'ConvergenceFactor'    : CONV - Alters the convergence criteria, 
%       M=M*CONV (default CONV = 1).
%    'DependentScattering'  : DEPSCAT - If false no near field dependency
%       for scattering between the cylinders will be regarded (default 
%       DEPSCAT = TRUE).
%
%   The amplitude scattering matrix T is calculated for a total of NANG 
%   scattering angles between 0 and 360 degrees and has dimension 2x2xNANG.
%   The considered angles are returned in vector ANG (in degrees).
%   
%   The cross sections (extinction, scattering and absorption) and 
%   wavenumber k are returned in the structure matrix C.
%
%
%   The calculation is based on the derivation of Lee [2]:
%   [2] Lee, S.-C., Dependent scattering of an obliquely incident plane 
%       wave by a collection of parallel cylinders. J. Appl. Phys. 68(10),
%       1990.
%
%   See also:
%   CALCCYL
%   CALCCYL_STRAT
%   CALCCYL_NF
%   CALCCYL_NF_MULTI
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Initialize parameters
p = inputParser;
p.addParamValue('HankelFunction', 1, @isscalar);
p.addParamValue('ConvergenceFactor', 1, @isscalar);
p.addParamValue('DependentScattering', true, @islogical);
p.parse(varargin{:});
        
hk = p.Results.HankelFunction;
conv = p.Results.ConvergenceFactor;
depscat = p.Results.DependentScattering;

if ~exist('iphi', 'var')
    iphi = 0;
end %if ~exist('iphi', 'var')
iphi = iphi/180*pi;     % the incident angle in radians

sgn = sign(hk-1.5);     % hk-convention dependent signum in exponentials

k = 2*pi/lambda*nm;     % wavenumber in outer medium nm
x = k*r;                % size parameter
m = ns/nm;              % relative refractive index

%% Calculate distance matrix
[xcg, ycg] = ndgrid(xc, yc);
dist = sqrt((xcg - xcg').^2 + (ycg - ycg').^2);

%% Calculate angle matrix
angm = acos((xcg - xcg')./dist)';
angm(isnan(angm)) = 0;
angm((ycg-ycg')<0) = -angm((ycg-ycg')<0);

%% Calculate phase shift vector
eps = exp(-sgn*1.i*k*(xc*cos(-iphi) + yc*sin(-iphi)));

%% Calculate expansion coefficients
[ann, bnp] = expcoeff_cyl_multi(x, m, k, dist, angm, eps, iphi, hk, ...
    conv, depscat);

%% Calculate amplitude scattering matrix
[T, C, ang] = asmcyl_multi(ann, bnp, k, dist, angm, nang, iphi, hk);
ang = ang/pi*180;

end
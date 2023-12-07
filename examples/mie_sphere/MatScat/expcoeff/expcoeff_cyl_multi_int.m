function [Ann, Bnp] = expcoeff_cyl_multi_int( ann, bnp, r, k, m, dist, ...
    angm, eps, iphi )
%EXPCOEFF_CYL_MULTI_INT Calculates the expansion coefficients of the 
%   internal field for the scattering by many infinite cylinders at
%   perpendicular incidence.
%
%   [Ann,Bnp] = EXPCOEFF_CYL_MULTI_INT(ann,bnp,R,K,M,DIST,ANGM,EPS,IPHI)
%   calculates the expansion coefficients Ann and Bnn for the internal 
%   fields of multiple cylinder scattering. The expansion coefficients of 
%   the external scattered field are submitted as ann and bnn. The 
%   cylinders have radius R, the wavenumber of the incoming light is K, the
%   relative refarctive index is NM. The cylinder positions are given in 
%   distance Matrix DIST and angle Matrix ANGM, EPS is the phase shift 
%   vector. The incident light is propagating in the xy-plane at an angle 
%   IPHI (in radians) according to the positive x-axis. The cylinder axis 
%   is aligned along the z-axis.
%
%   The calculation is based on the derivation of Lee [2]:
%   [2] Lee, S.-C., Dependent scattering of an obliquely incident plane 
%       wave by a collection of parallel cylinders. J. Appl. Phys. 68(10),
%       1990.
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Initialize parameters

hk = 2;             % allways use Hankel function of the 2nd kind
x = k*r;            % the size parameter
N = size(ann, 1);   % total number of cylidners

%% Get truncation number
MM = size(ann, 2);
M = (MM - 1)/2;
n = (-M:M)';
ng = meshgrid(n, 1:N);
s = (-M:M)';
[sg, ig] = meshgrid(n, 1:N);

%% Calculate the internal expansion coefficients
Bnp = eps(1:N).'*(besselj(n, x).'.*exp(1.j*n*iphi).') - ...
    bnp.*besselh(ng, hk, x);
Ann = eps(1:N).'*(besselj(n, x).'.*exp(1.j*n*iphi).') - ...
    ann.*besselh(ng, hk, x);

for i=1:N
    for n=-M:M
        [nu,z] = meshgrid((s-n)',k*dist(i,:)');
        bh = besselh(nu, hk, z);
        dum = bnp.*(-1.j).^(sg-n).*bh.*...
            exp(1.j*(s-n)*angm(:,i)').'.*besselj(n,x).*(ig~=i);
        Bnp(i,(n+M) + 1) = Bnp(i,(n+M) + 1) - sum(dum(isfinite(dum))); 
        dum = ann.*(-1.j).^(sg-n).*bh.*...
            exp(1.j*(s-n)*angm(:,i)').'.*besselj(n,x).*(ig~=i);
        Ann(i,(n+M) + 1) = Ann(i,(n+M) + 1) - sum(dum(isfinite(dum)));
    end %for n=-M:M
end %for i=1:N

Bnp = Bnp./besselj(ng, m*x)/m;
Ann = Ann./besselj(ng, m*x)/m^2;

end


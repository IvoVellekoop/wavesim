function [ann, bnp] = expcoeff_cyl_multi( x, m, k, dist, angm, eps, ...
    iphi, hk, conv, depscat )
%EXPCOEFF_CYL_MULTI Calculates the expansion coefficients for the
%   scattering by multiple infinite cylinders at perpendicular incidence.
%
%   [ANN,BNP]=EXPCOEFF_CYL_MULTI(X,M,K,DIST,ANGM,EPS,IPHI,HK,CONV,DEPSCAT) 
%   calculates the expansion coefficients ANN and BNP for multiple
%   cylinders with size parameters X and relative refractive indices m.
%   The wavenumber of the incident light is K. The cylinder positions are
%   given in distance Matrix DIST and angle Matrix ANGM, EPS is the phase 
%   shift vector. The incident light is propagating in the xy-plane at an 
%   angle IPHI (in radians) according to the positive x-axis. The cylinder 
%   axis is aligned along the z-axis. The Hankel function H_N^(HK) is used 
%   in the computation. The convergence criteria is multiplied by a factor 
%   of CONV. If DEPSCAT flag is set dependent scattering effects are
%   regarded.
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
if ~exist('hk', 'var')
    hk = 2;
end %if ~exist('hk', 'var')

if ~exist('conv', 'var')
    conv = 2;
end %if ~exist('conv', 'var')

if ~exist('depscat', 'var')
    depscat = true;
end %if ~exist('depscat', 'var')

N = numel(eps);         % total number of cylinders

sgn = sign(hk-1.5);     % hk-convention dependent signum in exponentials

%% Calculate truncation number
xmax = max(x);
M = ceil(conv*(xmax + 4*(xmax^(1/3)) + 2));

MM = 2*M + 1;
MN = N*MM;
n = -M:M;

%% Get single cylinder expansion coefficients
[~, an0, bn0, ~] = expcoeff_cyl(x, m, 90, hk, conv);
an0 = [xwrev(an0(2:end)) an0];
bn0 = [xwrev(bn0(2:end)) bn0];

%% Setup linear equation system
va = [];
vb = [];

for in=1:N
    va = [va an0*eps(in).*exp(sgn*1.i*n*iphi)];
    vb = [vb bn0*eps(in).*exp(sgn*1.i*n*iphi)];
end %for in=1:N

Sa = eye(MN);
Sb = eye(MN);

if depscat
    for l=1:N
        for s=-M:M
            for j=1:N
                for n=-M:M
                    if l~=j
                        Sa((j-1)*MM + (n+M) + 1, (l-1)*MM + (s+M) + 1) =...
                            (-sgn*1.i).^(s-n)*...
                            besselh(s-n, hk, k*dist(j,l))*...
                            exp(sgn*1.i*(s-n)*angm(l,j))*an0(n+M+1);
                        Sb((j-1)*MM + (n+M) + 1, (l-1)*MM + (s+M) + 1) =...
                            (-sgn*1.i).^(s-n)*...
                            besselh(s-n, hk, k*dist(j,l))*...
                            exp(sgn*1.i*(s-n)*angm(l,j))*bn0(n+M+1);
                    end %if l~=j
                end %for n=-M:M
            end %for j=1:N
        end %for s=-M:M
    end %for l=1:N
end %if depscat

%% Solve linear equation system
ann = linsolve(Sa,va.');
bnp = linsolve(Sb,vb.');

ann = reshape(ann, MM, N).';
bnp = reshape(bnp, MM, N).';
end


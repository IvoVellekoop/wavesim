function [T, C, ang] = asmcyl_multi( ann, bnp, k, dist, angm, nang, ...
    iphi, hk )
%ASMCYL_MULTI Calculates the amplitude scattering matrix and cross sections
%   for multiple cylindrical scatterers at perpendicular incidence.
%
%   [T,C,ANG] = ASMCYL_MULTI(ANN,BNP,K,DIST,ANGM,NANG,IPHI,HK) calculates 
%   the amplitude scattering matrix T and the cross sections C for the 
%   scattering of a monochromatic electromagnetic wave with wavenumber K by
%   many cylindrical particle with given expansion coefficients ANN and 
%   BNP. The cylinder positions are given in distance matrix DIST and angle
%   Matrix ANGM. The incident light is propagating in the xy-plane at an 
%   angle IPHI (in degree) according to the positive x-axis. Hankel 
%   function H_N^(HK) is used in the computation.
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
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Initialize parameters
sgn = sign(hk-1.5);     % hk-convention dependent signum in exponentials

%% Calculate truncation number
MM = size(ann,2);
M = (MM-1)/2;
n = -M:M;

%% Calculate amplitude scattering matrix
T = zeros(2,2,nang);
ang = (0:nang-1)/nang*2*pi;
for iang=1:nang
    expdum = exp(sgn*1.i*k*dist(1,:).*cos(angm(1,:) - ...
        (ang(iang) - iphi))).'*exp(sgn*1.i*n*(ang(iang)-iphi));
    T(1,1,iang) = sum(sum(bnp.*expdum));
    T(2,2,iang) = sum(sum(ann.*expdum));
end %for iang=1:nang

%% Calculate cross sections
dphi = ang(2)-ang(1);
C.ext(1) = 4/k*real(T(1,1,1));
C.ext(2) = 4/k*real(T(2,2,1));
C.sca(1) = 2/pi/k*sum(abs(T(1,1,:)).^2)*dphi;
C.sca(2) = 2/pi/k*sum(abs(T(2,2,:)).^2)*dphi;
C.abs(1) = C.ext(1) - C.sca(1);
C.abs(2) = C.ext(2) - C.sca(2);
C.k = k;


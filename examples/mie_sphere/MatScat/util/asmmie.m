function [S, C, ang] = asmmie( an, bn, nang, k )
%ASMMIE Calculates the amplitude scattering matrix and cross sections for 
%   given expansion coefficients of a spherical particle.
%
%   [S,C,ANG] = ASMMIE(AN,BN,NANG,K) calculates the amplitude scattering 
%   matrix S and the cross sections C for the scattering of a monochromatic
%   electromagnetic wave with wavenumber K by a spherical particle with 
%   given expansion coefficients AN and BN.
%
%   The amplitude scattering matrix S is calculated for a total of NANG 
%   scattering angles between 0 and 180 degrees and has dimension 2x2xNANG.
%   The considered angles are returned in vector ANG (in degrees).
%   
%   The cross sections (extinction, scattering and absorption) and 
%   wavenumber k are returned in the structure matrix C.
%
%   The calculation is based on the book of Bohren and Huffman [1]:
%   [1] Bohren, C. F. and Huffman, D. R., Absorption and scattering of 
%       light by small particles, Wiley-Interscience, New York, 1998.
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Initialize parameters
n = 1:numel(an);
ang = (0:nang-1)/(nang-1)*pi;

n2 = (2*n+1);
En = n2./(n.*(n+1));
anEn = an.*En;
bnEn = bn.*En;

%% Calculate theta dependent fucntions
[pin,taun] = angdepfun_mie(ang, n);

%% Calculate amplitude scattering matrix
S = zeros(2,2,nang);
S(1,1,:) = anEn*taun + bnEn*pin;
S(2,2,:) = anEn*pin + bnEn*taun;

%% Calculate cross sections
C.ext = 2*pi/k^2*sum(n2.*real(an + bn));
C.sca = 2*pi/k^2*sum(n2.*(abs(an).^2 + abs(bn).^2));
C.abs = C.ext - C.sca;
C.k = k;

end

function [T, C, ang] = asmcyl( anp, ann, bnp, bnn, nang, k )
%ASMCYL Calculates the amplitude scattering matrix and cross sections for 
%   given expansion coefficients of a single cylindrical particle.
%
%   [T,C,ANG] = ASMCYL(ANP,ANN,BNP,BNN,NANG,K) calculates the amplitude 
%   scattering matrix T and the cross sections C for the scattering of a 
%   monochromatic electromagnetic wave with wavenumber K by a cylindrical 
%   particle with given expansion coefficients ANP, ANN, BNP and BNN.
%
%   The amplitude scattering matrix T is calculated for a total of NANG 
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
n = 1:(numel(anp)-1);
ang = (0:nang-1)/(nang-1)*pi;

%% Calculate amplitude scattering matrix
T = zeros(2,2,nang);
for iang=1:nang
    T(1,1,iang) = bnp(1) + 2*sum(bnp(2:end).*cos(n*ang(iang)));
    T(2,2,iang) = ann(1) + 2*sum(ann(2:end).*cos(n*ang(iang)));
    T(2,1,iang) = -2*1.j*sum(anp(2:end).*sin(n*ang(iang)));
    T(1,2,iang) = -2*1.j*sum(bnn(2:end).*sin(n*ang(iang)));
end %for iang=1:nang

%% Calculate cross sections
C.ext(1) = 4/k*real(T(1,1,1));
C.ext(2) = 4/k*real(T(2,2,1));
C.sca(1) = 4/k*(abs(bnp(1)).^2 + 2*sum(abs(bnp(2:end)).^2));
C.sca(2) = 4/k*(abs(ann(1)).^2 + 2*sum(abs(ann(2:end)).^2));

C.abs(1) = C.ext(1) - C.sca(1);
C.abs(2) = C.ext(2) - C.sca(2);
C.k = k;

end

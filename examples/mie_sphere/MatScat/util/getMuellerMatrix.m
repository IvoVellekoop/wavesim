function [ M ] = getMuellerMatrix( S )
%GETMUELLERMATRIX Calculates the Mueller matrix for given amplitude
%   scattering matrix
%
%   [M] = GETMUELLERMATRIX(S) calculates the Mueller matrix elements
%   for given amplitude scattering matrix S. The resulting Matrix has a
%   size of 4x4 for each angle in S.
%
%   The calculation is based on the book of Bohren and Huffman [1]:
%   [1] Bohren, C. F. and Huffman, D. R., Absorption and scattering of 
%       light by small particles, Wiley-Interscience, New York, 1998.
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

M = zeros([4, 4, size(S(1,1,:))]);

S2 = S(1,1,:);
S4 = S(1,2,:);
S3 = S(2,1,:);
S1 = S(2,2,:);
S1abs = abs(S1).^2;
S2abs = abs(S2).^2;
S3abs = abs(S3).^2;
S4abs = abs(S4).^2;

M(1,1,:) = 0.5*(S1abs + S2abs + S3abs + S4abs);
M(1,2,:) = 0.5*(-S1abs + S2abs - S3abs + S4abs);
M(1,3,:) = real(S2.*conj(S3) + S1.*conj(S4));
M(1,4,:) = imag(S2.*conj(S3) - S1.*conj(S4));
M(2,1,:) = 0.5*(-S1abs + S2abs + S3abs - S4abs);
M(2,2,:) = 0.5*(S1abs + S2abs - S3abs - S4abs);
M(2,3,:) = real(S2.*conj(S3) - S1.*conj(S4));
M(2,4,:) = imag(S2.*conj(S3) + S1.*conj(S4));
M(3,1,:) = real(S2.*conj(S4) + S1.*conj(S3));
M(3,2,:) = real(S2.*conj(S4) - S1.*conj(S3));
M(3,3,:) = real(S1.*conj(S2) + S3.*conj(S4));
M(3,4,:) = imag(S2.*conj(S1) + S4.*conj(S3));
M(4,1,:) = imag(conj(S2).*S4 + conj(S3).*S1);
M(4,2,:) = imag(conj(S2).*S4 - conj(S3).*S1);
M(4,3,:) = imag(S1.*conj(S2) - S3.*conj(S4));
M(4,4,:) = real(S1.*conj(S2) - S3.*conj(S4));

end
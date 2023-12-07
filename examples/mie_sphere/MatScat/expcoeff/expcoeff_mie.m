function [an, bn] = expcoeff_mie( x, m, conv )
%EXPCOEFF_MIE Calculates the expansion coefficients for the Mie theory.
%
%   [AN,BN] = EXPCOEFF_MIE(X,M,CONV) calculates the expansion coefficients
%   AN and BN for a sphere with size parameter X and relative refractive 
%   index M. The convergence criteria is multiplied by a factor of CONV.
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
if ~exist('conv', 'var')
    conv = 1;
end %if ~exist('conv', 'var')

%% Calculate truncation number
M = ceil(conv*(x + 4*(x^(1/3)) + 2));
n = 1:M;

%% Calculate auxiliary parameters
Sx = ricbesj(n, x);
dSx = dricbesj(n, x);
Smx = ricbesj(n, m*x);
dSmx = dricbesj(n, m*x);
xix = ricbesh(n, x);
dxix = dricbesh(n, x);

%% Calculate expansion coefficients
an = (m*Smx.*dSx - Sx.*dSmx) ./ (m*Smx.*dxix - xix.*dSmx);
bn = (Smx.*dSx - m*Sx.*dSmx) ./ (Smx.*dxix - m*xix.*dSmx);

end

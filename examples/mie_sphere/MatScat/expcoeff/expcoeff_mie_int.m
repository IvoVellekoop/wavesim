function [fn, gn] = expcoeff_mie_int( x, m, N )
%EXPCOEFF_MIE Calculates the expansion coefficients of the internal fields
%   for the Mie theory.
%
%   [FN,GN] = EXPCOEFF_MIE(X,M,N) calculates the internal fields 
%   expansion coefficients FN and GN for a sphere with size parameter X and
%   relative refractive index M. A total of N coefficients are given.
%
%   The calculation is based on the book of Bohren and Huffman [1]:
%   [1] Bohren, C. F. and Huffman, D. R., Absorption and scattering of 
%       light by small particles, Wiley-Interscience, New York, 1998.
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Calculate truncation number
n = (1:N)';

%% Calculate auxiliary parameters
Sx = ricbesj(n, x);
dSx = dricbesj(n, x);
Smx = ricbesj(n, m*x);
dSmx = dricbesj(n, m*x);
xix = ricbesh(n, x);
dxix = dricbesh(n, x);

%% Calculate expansion coefficients
fn = (m*Sx.*dxix - m*xix.*dSx) ./ (Smx.*dxix - m*xix.*dSmx);
gn = (m*Sx.*dxix - m*xix.*dSx) ./ (m*Smx.*dxix - xix.*dSmx);

end

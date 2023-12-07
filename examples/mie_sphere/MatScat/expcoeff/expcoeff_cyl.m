function [anp, ann, bnp, bnn] = expcoeff_cyl( x, m, zeta, hk, conv )
%EXPCOEFF_CYL Calculates the expansion coefficients for the infinite 
%   circular cylinder.
%
%   [ANP,ANN,BNP,BNN] = EXPCOEFF_CYL(X,M,ZETA,HK,CONV) calculates the 
%   expansion coefficients ANP,ANN,BNP and BNN for a cylinder with size 
%   parameter X and relative refractive index M. The cylinder axis is 
%   aligned along the z-axis and rotated by an angle ZETA (in degree)
%   around the y-axis. The incident light is propagating along the positive
%   z-axis. If not defined ZETA is set to 90 degrees (perpendicular 
%   incidence). The Hankel function H_N^(HK) is used in the computation. 
%   The convergence criteria is multiplied by a factor of CONV.
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
if ~exist('zeta', 'var')
    zeta = 90;
end %if ~exist('zeta', 'var')

zeta = zeta/180*pi;     % inclination angle in radians

if ~exist('hk', 'var')
    hk = 1;
end %if ~exist('hk', 'var')

if ~exist('conv', 'var')
    conv = 1;
end %if ~exist('conv', 'var')

%% Calculate truncation number
xsin = x*sin(zeta);
M = ceil(conv*(xsin + 4*(xsin^(1/3)) + 2));
n = 0:M;

%% Calculate auxiliary variables
xi = x*sin(zeta);
eta = x*sqrt(m^2 - cos(zeta)^2);
jneta = besselj(n, eta);
djneta = dbesselj(n, eta);
jnxi = besselj(n, xi);
djnxi = dbesselj(n, xi);
hnxi = besselh(n, hk, xi);
dhnxi = dbesselh(n, hk, xi);

An = 1.j*xi*(xi*djneta.*jnxi - eta*jneta.*djnxi);
Bn = xi*(m^2*xi*djneta.*jnxi - eta*jneta.*djnxi);
Cn = n*cos(zeta)*eta.*jneta.*jnxi*(xi^2/eta^2 - 1);
Dn = n*cos(zeta)*eta.*jneta.*hnxi*(xi^2/eta^2 - 1);
Vn = xi*(m^2*xi*djneta.*hnxi - eta*jneta.*dhnxi);
Wn = 1.j*xi*(eta*jneta.*dhnxi - xi*djneta.*hnxi);

wvd = (Wn.*Vn + 1.j*Dn.^2);
cd = 1.j*Cn.*Dn;

idx = isnan(wvd);

%% Calculate expansion coefficients
anp = (Cn.*Vn - Bn.*Dn)./wvd;
ann = -(An.*Vn - cd)./wvd;
bnp = (Wn.*Bn + cd)./wvd;
bnn = -1.j*(Cn.*Wn + An.*Dn)./wvd;

anp(idx) = 0.;
ann(idx) = 0.;
bnp(idx) = 0.;
bnn(idx) = 0.;
end

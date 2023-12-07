function [ E, H, P, S, C, ang ] = calcmie_nf( r, ns, nm, lambda, ...
    xc, yc, zc, varargin )
%CALCMIE_NF Calculate the near field solution for the scattering of 
%   electromagnetic radiation by a single (stratified) sphere.
%
%   [E,H,P,S,C,ANG] = CALCMIE_NF(R,NS,NM,LAMBDA,XC,YC,ZC,VARARGIN) 
%   calculates the near E- and H-fields for the scattering of a 
%   monochromatic electromagnetic wave with wavelength LAMBDA by a sphere 
%   with radius R and (complex) refractive index NS. NM is the refractive 
%   index of the surrounding medium. The incident light is propagating 
%   along the z-direction. The E-field of the incoming electromagnetic wave
%   is polarized along the x-direction. The near field solution is 
%   calculated for 3D cartesian coordinates given in arrays XC, YC and ZC.
%   The far field solution is returned in [S,C,ANG].
%
%   For a stratified sphere R and NS are given as one-dimensional lists. 
%   The radii and refractive index lists have to be the same size, the 
%   radii have to be sorted from smallest to largest and the corresponding
%   refractive indices have to be in the same order.
%
%   VARARGIN can take the following keywords:
%    'ConvergenceFactor': CONV - Alters the convergence criteria M=M*CONV
%                         (default CONV = 1).
%    'TotalField'       : TF_FLAG - If set calculate total fields, 
%                         otherwise scattered fields are returned (default
%                         TF_FLAG = FALSE).
%    'Cartesian'        : CC_FLAG - get the result in cartesian coordinates
%                         (default CC_FLAG = True).
%    'nang'             : NANG - number of far field angles to evaluate 
%                         (default NANG = 1800).
%
%   The calculation is based on the book of Bohren and Huffman [1]:
%   [1] Bohren, C. F. and Huffman, D. R., Absorption and scattering of 
%       light by small particles, Wiley-Interscience, New York, 1998.
%
%   SYNTAX:
%
%   [E, H] = calcmie_nf( r, ns, nm, lambda, xc, yc, zc )
%   [E, H, P] = calcmie_nf( r, ns, nm, lambda, xc, yc, zc )
%   [E, H, ~, S, C, ang] = calcmie_nf( r, ns, nm, lambda, xc, yc, zc )
%
%   See also:
%   CALCMIE
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Initialize parameters
p = inputParser;
p.addParamValue('ConvergenceFactor', 1, @isscalar);
p.addParamValue('TotalField', false, @islogical);
p.addParamValue('Cartesian', true, @islogical);
p.addParamValue('nang', 1800, @isscalar);
p.parse(varargin{:});
    
conv = p.Results.ConvergenceFactor;
tf_flag = p.Results.TotalField;
cc_flag = p.Results.Cartesian;

nang = p.Results.nang;

strtfd = false;
if numel(r) > 1
    strtfd = true;
end %if numel(r) > 1

if numel(r) ~= numel(ns)
    error('Radius list r and refractive index list NS have to be same size!');
end %if numel(r) ~= numel(n)

k = 2*pi/lambda*nm; % the wavenumber in medium nm
x = k*r;            % the size parameter
m = ns/nm;          % the relative refractive index

%% Calculate expansion coefficients
if strtfd
    [an, bn] = expcoeff_mie_strat(x, m, conv);
else %strtfd
    [an, bn] = expcoeff_mie(x, m, conv);
end %strtfd

%% Calculate field solution
[E,H] = nfmie(an, bn, xc, yc, zc, r, ns, nm, lambda, tf_flag, cc_flag);

%% Calculate Poynting vector
if nargout > 2
    P = real(cross(E, conj(H)));
end %nargout > 2

%% Calculate Far Field solution
if nargout > 3    
    [S, C, ang] = asmmie(an, bn, nang, k);
    ang = ang/pi*180;
end %nargout > 3

end

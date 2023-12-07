function [T, C, ang] = calccyl( r, ns, nm, lambda, nang, zeta, varargin )
%CALCCYL Calculates the amplitude scattering matrix and cross sections for
%   the scattering of electromagnetic radiation by a single infinite 
%   (stratified) cylinder at oblique incidence.
%
%   [T,C,ANG] = CALCCYL(R,NS,NM,LAMBDA,NANG,ZETA,VARARGIN) calculates the 
%   amplitude scattering matrix T and the cross sections C for the 
%   scattering of a monochromatic electromagnetic wave with vacuum 
%   wavelength LAMBDA by an infinite circular cylinder with radius R and 
%   (complex) refractive index NS. NM is the refractive index of the 
%   surrounding medium. The cylinder axis is aligned along the z-axis and 
%   rotated by an angle ZETA (in degree) around the y-axis. The incident 
%   light is propagating along the positive z-axis. If not defined ZETA is 
%   set to 90 degrees (perpendicular incidence).
%
%   For a stratified cylinder R and NS are given as one-dimensional lists. 
%   The radii and refractive index lists have to be the same size, the 
%   radii have to be sorted from smallest to largest and the corresponding
%   refractive indices have to be in the corresponding order.
%
%   VARARGIN can take the following keywords:
%    'ConvergenceFactor'    : CONV - Alters the convergence criteria, 
%       M=M*CONV (default CONV = 1).
%
%   The amplitude scattering matrix T is calculated for a total of NANG 
%   scattering angles between 0 and 180 degrees and has dimension 2x2xNANG.
%   The considered angles are returned in vector ANG (in degrees).
%   
%   The cross sections (extinction, scattering and absorption) and 
%   wavenumber k are returned in the structure matrix C.
%
%   The basic calculation is based on the book of Bohren and Huffman [1],
%   the calculation for the stratified cylinder is based on the book of 
%   Kerker [5].
%   [1] Bohren, C. F. and Huffman, D. R., Absorption and scattering of 
%       light by small particles, Wiley-Interscience, New York, 1998.
%   [5] Kerker, M., The scattering of light and other electromagnetic 
%       radiation, Academic Press, 1969
%
%   SYNTAX:
%
%   [T, C, ANG] = calccyl( r, ns, nm, lambda, nang )
%   [T, C, ANG] = calccyl( r, ns, nm, lambda, nang, zeta )
%
%   See also:
%   CALCCYL_MULTI
%   CALCCYL_NF
%   CALCCYL_NF_MULTI
%
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Initialize parameters
p = inputParser;
p.addParamValue('ConvergenceFactor', 1, @isscalar);
p.parse(varargin{:});
        
conv = p.Results.ConvergenceFactor;

strtfd = false;
if numel(r) > 1
    strtfd = true;
end %if numel(r) > 1

if numel(r) ~= numel(ns)
    error('Radius list and refractive index list have to be same size!');
end %if numel(r) ~= numel(n)

if ~exist('zeta', 'var')
    zeta = 90;
end %if ~exist('zeta', 'var')

if strtfd && zeta ~= 90
    error(['Stratified cylinder solution is only applicable for', ...
        'perpendicular incidence (zeta = 90 degrees)!']);
end %if strtfd && zeta ~= 90

k = 2*pi/lambda*nm;     % wavenumber in outer medium nm
x = k*r;                % size parameter
m = ns/nm;              % relative refractive index

%% Calculate expansion coefficients
if strtfd
    [ann, bnp] = expcoeff_cyl_strat(x, m, 1, conv);
    anp = zeros(size(ann));
    bnn = zeros(size(ann));
else %strtfd
    [anp, ann, bnp, bnn] = expcoeff_cyl(x, m, zeta, 1, conv);
end %strtfd

%% Calculate amplitude scattering matrix and cross sections
[T, C, ang] = asmcyl(anp, ann, bnp, bnn, nang, k);
ang = ang/pi*180;

end
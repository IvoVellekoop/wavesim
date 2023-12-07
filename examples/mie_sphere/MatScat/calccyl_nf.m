function [ Ep, En, Hp, Hn, Sp, Sn, T, C, ang ] = calccyl_nf( r, ns, nm, ...
    lambda, xc, yc, zc, zeta, varargin )
%CALCCYL_NF Calculate the near field solution for the scattering of 
%   electromagnetic radiation by an infinite (stratified) cylinder.
%
%   [EP,EN,HP,HN,SP,SN,T,C,ANG] = CALCCYL_NF(R,NS,NM,LAMBDA,XC,YC,ZC,ZETA,
%   VARARGIN) calculates the near E- and H-fields for the scattering of a 
%   monochromatic electromagnetic wave with vacuum wavelength LAMBDA by a 
%   (stratified) cylinder with radius R and (complex) refractive index NS. 
%   NM is the refractive index of the surrounding medium. The near field 
%   solution is calculated for 3D cartesian coordinates given in arrays 
%   XC, YC and ZC. The cylinder axis is aligned along the z-axis and 
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
%   The field solution is for incoming electromagnetic waves polarized 
%   perpendicular (EP,HP) and normal (EN,HN) according to the cylinder 
%   axis. Additionally the Poynting vectors are returned for both 
%   polarizations (SP,SN). The dimension of the fields and Poynting vectors
%   is size(XC)x3. The far field solution is returned in [T,C,ANG].
%
%   The calculation is based on the book of Bohren and Huffman [1]:
%   [1] Bohren, C. F. and Huffman, D. R., Absorption and scattering of 
%       light by small particles, Wiley-Interscience, New York, 1998.
%
%
%   SYNTAX:
%
%   [Ep, En, Hp, Hn] = calccyl_nf(r, ns, nm, lambda, xc, yc, zc);
%   [Ep, En, Hp, Hn, Sp, Sn] = calccyl_nf(r, ns, nm, lambda, ...
%                              xc, yc, zc, zeta);
%   [Ep, En, Hp, Hn, Sp, Sn, T, C, ang] = calccyl_nf(r, ns, nm, lambda, ...
%                                         xc, yc, zc, zeta, 'nang', nang);
%
%
%   See also:
%   CALCCYL
%   CALCCYL_MULTI
%   CALCCYL_NF_MULTI
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

hk = 1; % allways use Hankel function of the 1st kind

strtfd = false;
if numel(r) > 1
    strtfd = true;
end %if numel(r) > 1

if ~exist('zeta', 'var')
    zeta = 90;
end %if ~exist('zeta', 'var')

if strtfd && zeta ~= 90
    error(['Stratified cylinder solution is only applicable for ', ...
        'perpendicular incidence (zeta = 90 degrees)!']);
end %if strtfd && zeta ~= 90

k = 2*pi/lambda*nm; % the wavenumber in medium nm
x = k*r;            % the size parameter
m = ns/nm;          % the relative refractive index

%% Calculate expansion coefficients
if strtfd
    [ann, bnp] = expcoeff_cyl_strat(x, m, 1, conv);
    anp = zeros(size(ann));
    bnn = zeros(size(ann));
else %strtfd
    [anp, ann, bnp, bnn] = expcoeff_cyl(x, m, zeta, hk, conv);
end %strtfd

%% Calculate field solution
[Ep, En, Hp, Hn] = nfcyl(anp, ann, bnp, bnn, xc, yc, zc, ...
    r, ns, nm, lambda, zeta, tf_flag, cc_flag);

%% Calculate Poynting vector
if nargout > 4
    Sp = real(cross(Ep, conj(Hp)));
    Sn = real(cross(En, conj(Hn)));
end %nargout > 4

%% Calculate Far Field solution
if nargout > 6    
    [T, C, ang] = asmcyl(anp, ann, bnp, bnn, nang, k);
    ang = ang/pi*180;
end %nargout > 6

end

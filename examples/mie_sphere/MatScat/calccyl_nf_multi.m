function [ Ep, En, Hp, Hn, Sp, Sn, T, C, ang ] = calccyl_nf_multi( ...
    r, ns, nm, lambda, xc, yc, xf, yf, varargin )
%CALCCYL_NF_MULTI Calculate the near field solution for the scattering of 
%   electromagnetic radiation by multiple infinite cylinders at 
%   perpendicular incidence.
%
%   [EP,EN,HP,HN,SP,SN,T,C,ANG] = CALCCYL_NF_MULTI(R,NS,NM,LAMBDA,XC,YC,XF,
%   YF,VARARGIN) calculates the near E- and H-fields for the scattering of 
%   a monochromatic electromagnetic wave with vacuum wavelength LAMBDA at 
%   perpendicular incidence by many equal cylinders with radius R and 
%   (complex) refractive index NS. NM is the refractive index of the 
%   surrounding medium. The coordinates of the cylinders are given in XC 
%   and YC, the cylinders may not overlap. The near field solution is 
%   calculated for 2D cartesian coordinates given in arrays XF and YF. The 
%   coordinate fields XC and YC as well as XF and YF have to be the same 
%   size.
%
%   VARARGIN can take the following keywords:
%    'ConvergenceFactor': CONV - Alters the convergence criteria M=M*CONV
%                         (default CONV = 1).
%    'TotalField'       : TF_FLAG - If set calculate total fields, 
%                         otherwise scattered fields are returned (default
%                         TF_FLAG = FALSE).
%    'PoyntingDirect'   : PD_FLAG - calculate the Poynting vector using the
%                         analytical formula instead of performing the
%                         cross product of the calculated fields (default
%                         PD_FLAG = False).
%    'Cartesian'        : CC_FLAG - get the result in cartesian coordinates
%                         (default CC_FLAG = True).
%    'nang'             : NANG - number of far field angles to evaluate 
%                         (default NANG = 1800).
%    'iphi'             : IPHI - incident phi angle in degrees (default
%                         IPHI = 0).
%
%   The field solution is for incoming electromagnetic waves polarized 
%   perpendicular (EP,HP) and normal (EN,HN) according to the cylinder 
%   axis. Additionally the Poynting vectors are returned for both 
%   polarizations (SP,SN). The dimension of the fields and Poynting vectors
%   is size(XC)x3. The far field solution is returned in [T,C,ANG].
%
%
%   The calculation is based on the derivation of Lee [2]:
%   [2] Lee, S.-C., Dependent scattering of an obliquely incident plane 
%       wave by a collection of parallel cylinders. J. Appl. Phys. 68(10),
%       1990.
%
%
%   SYNTAX:
%
%   [Ep, En, Hp, Hn] = calccyl_nf_multi(r, ns, nm, lambda, xc, yc, xf, yf);
%   [Ep, En, Hp, Hn, Sp, Sn] = calccyl_nf_multi(r, ns, nm, lambda, ...
%                              xc, yc, xf, yf);
%   [Ep, En, Hp, Hn, Sp, Sn, T, C, ang] = calccyl_nf_multi(r, ns, nm, ...
%                                         lambda, xc, yc, xf, yf, ...
%                                         'nang', nang);
%
%
%   See also:
%   CALCCYL
%   CALCCYL_STRAT
%   CALCCYL_MULTI
%   CALCCYL_NF
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Initialize parameters
p = inputParser;
p.addParamValue('ConvergenceFactor', 1, @isscalar);
p.addParamValue('TotalField', false, @islogical);
p.addParamValue('PoyntingDirect', false, @islogical);
p.addParamValue('Cartesian', true, @islogical);
p.addParamValue('nang', 1800, @isscalar);
p.addParamValue('iphi', 0, @isscalar);
p.parse(varargin{:});

conv = p.Results.ConvergenceFactor;
tf_flag = p.Results.TotalField;
pd_flag = p.Results.PoyntingDirect;
cc_flag = p.Results.Cartesian;
nang = p.Results.nang;
iphi = p.Results.iphi/180*pi;

hk = 2; % allways use Hankel function of the 2nd kind

% Definition of refractive index is: n = nr - i*ni
% WARNING: for solutions using Hankel funtion of the 2nd kind
if imag(ns) > 0.
    ns = real(ns) - imag(ns)*1.j;
end %if imag(n_s) > 1

k = 2*pi/lambda*nm; % the wavenumber in medium nm
x = k*r;            % the size parameter
m = ns/nm;          % the relative refractive index

%% Calculate distance matrix
[xcg, ycg] = ndgrid(xc, yc);
dist = sqrt((xcg - xcg').^2 + (ycg - ycg').^2);

%% Calculate angle matrix
ang = acos((xcg - xcg')./dist)';
ang(isnan(ang)) = 0;
ang((ycg - ycg')<0) = -ang((ycg - ycg')<0);

%% Calculate phase shift vector
eps = exp(-1.i*k*(xc*cos(iphi) - yc*sin(iphi)));

%% Calculate expansion coefficients

% external fields
[ann, bnp] = expcoeff_cyl_multi(x, m, k, dist, ang, eps, iphi, hk, conv);

% internal fields
[Ann, Bnp] = expcoeff_cyl_multi_int(ann, bnp, r, k, m, dist, ang, eps, iphi);

%% Calculate field solution
[Ep, En, Hp, Hn] = nfcyl_multi(ann, bnp, Ann, Bnp, r, ns, nm, lambda, ...
    xc, yc, xf, yf, iphi, tf_flag, cc_flag);

%% Calculate Poynting vector
if nargout > 4
    if pd_flag
        [Sp, Sn] = nfcyl_multi_Poynting(ann, bnp, Ann, Bnp, r, ns, nm, ...
            lambda, xc, yc, xf, yf, tf_flag, cc_flag);        
    else %if pd_flag
        Sp = real(cross(Ep, conj(Hp)));
        Sn = real(cross(En, conj(Hn)));
    end %if pd_flag
end %nargout > 4


%% Calculate Far Field solution
if nargout > 6    
    [T, C, ang] = asmcyl_multi(ann, bnp, k, dist, ang, nang, iphi, hk);
    ang = ang/pi*180;
end %nargout > 6

end


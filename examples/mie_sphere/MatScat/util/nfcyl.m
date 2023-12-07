function [ Ep, En, Hp, Hn ] = nfcyl( anp, ann, bnp, bnn, xc, yc, zc, ...
    r, ns, nm, lambda, zeta, tf_flag, cc_flag )
%NFCYL Calculates the near field solution for given expansion
%   coefficients of a cylindrical particle.
%
%   [EP,EN,HP,HN] = NFCYL(ANP,ANN,BNP,BNN,XC,YC,ZC,R,NS,NM,LAMBDA,ZETA,
%   TF_FLAG, CC_FLAG) calculates the near field E and H for the scattering
%   of a monochromatic electromagnetic wave by a cylindrical particle with
%   given expansion coefficients ANP, ANN, BNP and BN. The solution is
%   for incoming electromagnetic waves polarized perpendicular (EP,HP)
%   and normal (EN, HN) according to the cylinder axes. The near field
%   solution is calculated for 3D cartesian coordinates given in arrays
%   XC, YC and ZC. The outer radius of the particle is given by R, the
%   (complex) refractive index of the cylinder is NS. LAMBDA is the
%   wavelength of the incident light in vacuum, NM is the refractive index
%   of the outer medium. The cylinder axis is aligned along the z-axis and
%   rotated by an angle ZETA (in degree) around the y-axis. The incident
%   light is propagating along the positive z-axis. If not defined ZETA is
%   set to 90 degrees (perpendicular incidence).
%
%   For a stratified cylinder R and NS are given as one-dimensional lists. 
%   The radii and refractive index lists have to be the same size, the 
%   radii have to be sorted from smallest to largest and the corresponding
%   refractive indices have to be in the corresponding order.
%
%   If TF_FLAG is set total fields are returned, otherwise scattered fields
%   are calculated. IF CC_FLAG is set resulting vectors are in cartesian
%   coordinate vector basis, in cylindrical coordinate basis otherwise.
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

if ~exist('tf_flag', 'var')
    tf_flag = false;
end %if ~exist('tf_flag', 'var')

if ~exist('cc_flag', 'var')
    cc_flag = false;
end %if ~exist('cc_flag', 'var')

MatScat_const;

strtfd = false;
if numel(r) > 1
    strtfd = true;
end %if numel(r) > 1

if strtfd && zeta ~= 90
    error(['Stratified cylinder solution is only applicable for', ...
        'perpendicular incidence (zeta = 90 degrees)!']);
end %if strtfd && zeta ~= 90

if zeta ~= 90
    warning(['Calculation of internal fields is only applicable for ', ...
        'perpendicular incidence (zeta = 90 degrees)!']);
end %if strtfd && zeta ~= 90

zeta = zeta/180*pi;         % inclination angle in radians

k = 2*pi/lambda*nm;         % the wavenumber in medium nm
x = k*r;                    % the size parameter
m = ns/nm;                  % the relative refractive index

k1 = 2*pi/lambda*ns;        % the wavenumber in medium ns
omega = 2.*pi/lambda*c0;    % the angular frequency

%% Get truncation number
M = numel(anp)-1;
n = (-M:M)';

%% Get scattered field expansion coefficients
anp = [-xwrev(anp(2:end)) anp].';
ann = [xwrev(ann(2:end)) ann].';
bnp = [xwrev(bnp(2:end)) bnp].';
bnn = [-xwrev(bnn(2:end)) bnn].';

%% Get internal field expansion coefficients
if strtfd   
    [fnp, gnn, vnp, wnn] = expcoeff_cyl_strat_int( ann, bnp, r, k, m );

    fnn = zeros(size(fnp));
    gnp = zeros(size(fnp));
    vnn = zeros(size(fnp));
    wnp = zeros(size(fnp));
else %if strtfd
    % TODO: oblique angle
    fnp = (besselj(n, x) - bnp.*besselh(n, x))./besselj(n, m*x)/m;
    gnn = (besselj(n, x) - ann.*besselh(n, x))./besselj(n, m*x)/m^2;
    fnn = zeros(size(fnp));
    gnp = zeros(size(fnp));
end %if strtfd

%% Calculate incident field expansion coefficient
E0 = 1;
E_n = E0*(-1.j).^n/k/sin(zeta);

%% Define field arrays
Ep = zeros([numel(xc), 3]);
En = zeros([numel(xc), 3]);
Hp = zeros([numel(xc), 3]);
Hn = zeros([numel(xc), 3]);

for ic=1:numel(xc)
    
    % get polar coordinates
    [phi, rad, z] = cart2pol(xc(ic), yc(ic), zc(ic));
    
    % calculate coordinate transformation matrix
    R = getTransformationMatrix(phi);
    
    % Incident Field propagating along e_x (not -e_x) for zeta=90°
    phi = phi + pi;
    
    % define arrays
    Mn = zeros(2*M+1,3);
    Nn = zeros(2*M+1,3);
    
    Edump = zeros(2*M+1,3);
    Edumn = zeros(2*M+1,3);
    Hdump = zeros(2*M+1,3);
    Hdumn = zeros(2*M+1,3);
    
    if rad > r(end)
        % scattered field
        
        % auxiliary variables
        h = -k*cos(zeta);
        skh = sqrt(k^2 - h^2);
        rho = rad*skh;
        hn = besselh(n,rho);
        dhn = dbesselh(n,rho);
        
        Mn(:,1) = 1.j*n/rho.*hn;
        Mn(:,2) = -dhn;
        Nn(:,1) = 1.j*h*dhn;
        Nn(:,2) = -h/rho*n.*hn;
        Nn(:,3) = skh*hn;
        for d=1:3
            Mn(:,d) = Mn(:,d)*skh.*exp(1.j*(n.*phi + h*z));
            Nn(:,d) = Nn(:,d)*skh/k.*exp(1.j*(n.*phi + h*z));
            
            Edump(:,d) = -E_n.*(bnp.*Nn(:,d) + 1.j*anp.*Mn(:,d));
            Hdump(:,d) = E_n.*(bnp.*Mn(:,d) + 1.j*anp.*Nn(:,d));
            Edumn(:,d) = E_n.*(1.j*ann.*Mn(:,d) + bnn.*Nn(:,d));
            Hdumn(:,d) = -E_n.*(1.j*ann.*Nn(:,d) + bnn.*Mn(:,d));
        end %for d=1:3
        Ep(ic,:) = sum(Edump);
        En(ic,:) = sum(Edumn);
        Hp(ic,:) = 1.j*k/omega/mue0*sum(Hdump);
        Hn(ic,:) = 1.j*k/omega/mue0*sum(Hdumn);
        
        Ep(ic,:) = R*Ep(ic,:).';
        En(ic,:) = R*En(ic,:).';
        Hp(ic,:) = R*Hp(ic,:).';
        Hn(ic,:) = R*Hn(ic,:).';
        
        if tf_flag
            % total field
            Ep(ic,3) = Ep(ic,3) + exp(1.j*k*xc(ic));
            Hp(ic,2) = Hp(ic,2) - nm/c0/mue0*exp(1.j*k*xc(ic));
            En(ic,2) = En(ic,2) - exp(1.j*k*xc(ic));
            Hn(ic,3) = Hn(ic,3) - nm/c0/mue0*exp(1.j*k*xc(ic));
        end %if tf_flag
        
        if ~cc_flag
            Ep(ic,:) = R'*Ep(ic,:).';
            En(ic,:) = R'*En(ic,:).';
            Hp(ic,:) = R'*Hp(ic,:).';
            Hn(ic,:) = R'*Hn(ic,:).';
        end %if ~cc_flag
        
    else %if rad > r
        ilay = max([find(r >=rad, 1, 'first') , find(r < rad, 1, 'last')]);
        % internal total field

        % auxiliary variables
        h = -k1(ilay)*cos(zeta);
        skh = sqrt(k1(ilay)^2 - h^2);
        rho = rad*skh;
        zn = besselj(n,rho);
        dzn = dbesselj(n,rho);
        
        Mn(:,1) = 1.j*n/rho.*zn;
        Mn(:,2) = -dzn;
        Nn(:,1) = 1.j*h*dzn;
        Nn(:,2) = -h/rho*n.*zn;
        Nn(:,3) = skh*zn;
        for d=1:3
            Mn(:,d) = Mn(:,d)*skh.*exp(1.j*(n.*phi + h*z));
            Nn(:,d) = Nn(:,d)*skh/k1(ilay).*exp(1.j*(n.*phi + h*z));
            
            Edump(:,d) = E_n.*(fnp(:,ilay).*Nn(:,d) + ...
                gnp(:,ilay).*Mn(:,d));
            Hdump(:,d) = -E_n.*(fnp(:,ilay).*Mn(:,d) + ...
                gnp(:,ilay).*Nn(:,d));
            Edumn(:,d) = -E_n.*(gnn(:,ilay).*Mn(:,d) + ...
                fnn(:,ilay).*Nn(:,d));
            Hdumn(:,d) = -E_n.*(gnn(:,ilay).*Nn(:,d) + ...
                fnn(:,ilay).*Mn(:,d));
        end
        
        if strtfd
            zn = besselh(n,rho);
            dzn = dbesselh(n,rho);
            
            Mn(:,1) = 1.j*n/rho.*zn;
            Mn(:,2) = -dzn;
            Nn(:,1) = 1.j*h*dzn;
            Nn(:,2) = -h/rho*n.*zn;
            Nn(:,3) = skh*zn;
            for d=1:3
                Mn(:,d) = Mn(:,d)*skh.*exp(1.j*(n.*phi + h*z));
                Nn(:,d) = Nn(:,d)*skh/k1(ilay).*exp(1.j*(n.*phi + h*z));
                
                Edump(:,d) = Edump(:,d) - ...
                    E_n.*(vnp(:,ilay).*Nn(:,d) + wnp(:,ilay).*Mn(:,d));
                Hdump(:,d) = Hdump(:,d) + ...
                    E_n.*(vnp(:,ilay).*Mn(:,d) + wnp(:,ilay).*Nn(:,d));
                Edumn(:,d) = Edumn(:,d) + ...
                    E_n.*(wnn(:,ilay).*Mn(:,d) + vnn(:,ilay).*Nn(:,d));
                Hdumn(:,d) = Hdumn(:,d) + ...
                    E_n.*(wnn(:,ilay).*Nn(:,d) + vnn(:,ilay).*Mn(:,d));
            end
        end %if strtfd
        
        Ep(ic,:) = sum(Edump);
        En(ic,:) = 1.j*sum(Edumn);
        Hp(ic,:) = 1.j*k1(ilay)/omega/mue0*sum(Hdump);
        Hn(ic,:) = k1(ilay)/omega/mue0*sum(Hdumn);
        
        Ep(ic,:) = R*Ep(ic,:).';
        En(ic,:) = R*En(ic,:).';
        Hp(ic,:) = R*Hp(ic,:).';
        Hn(ic,:) = R*Hn(ic,:).';
        
        if ~tf_flag
            % internal scattered field
            Ep(ic,3) = Ep(ic,3) - exp(1.j*k*xc(ic));
            Hp(ic,2) = Hp(ic,2) + nm/c0/mue0*exp(1.j*k*xc(ic));
            En(ic,2) = En(ic,2) + exp(1.j*k*xc(ic));
            Hn(ic,3) = Hn(ic,3) + nm/c0/mue0*exp(1.j*k*xc(ic));
        end %if ~tf_flag
        
        if ~cc_flag
            Ep(ic,:) = R'*Ep(ic,:).';
            En(ic,:) = R'*En(ic,:).';
            Hp(ic,:) = R'*Hp(ic,:).';
            Hn(ic,:) = R'*Hn(ic,:).';
        end %if ~cc_flag
               
    end %if rad > r
    
end %for icoord=1:length(coords)

Ep = reshape(Ep,[size(xc), 3]);
En = reshape(En,[size(xc), 3]);
Hp = reshape(Hp,[size(xc), 3]);
Hn = reshape(Hn,[size(xc), 3]);

end

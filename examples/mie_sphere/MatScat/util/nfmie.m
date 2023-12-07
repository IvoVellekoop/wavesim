function [ E, H ] = nfmie( an, bn, xf, yf, zf, r, ns, nm, lambda, ...
    tf_flag, cc_flag )
%NFMIE Calculates the near field solution for given expansion
%   coefficients.
%
%   [E,H] = NFMIE(AN,BN,XC,YC,ZC,R,NM,LAMBDA,TF_FLAG,CC_FLAG) calculates
%   the near field E and H for the scattering of a monochromatic
%   electromagnetic wave by a particle with given expansion coefficients
%   AN and BN. The near field solution is calculated for 3D cartesian
%   coordinates given in arrays XC, YC and ZC. The outer radius of the
%   particle is given by R, LAMBDA is the wavelength of the incident light,
%   NM is the refractive index of the outer medium.
%
%   For a stratified sphere R and NS are given as one-dimensional lists. 
%   The radii and refractive index lists have to be the same size, the 
%   radii have to be sorted from smallest to largest and the corresponding
%   refractive indices have to be in the same order.
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

k = 2*pi/lambda*nm;         % the wavenumber in medium nm
x = k*r;                    % the size parameter
m = ns/nm;                  % the relative refractive index

k1 = 2*pi/lambda*ns;        % the wavenumber in medium ns
omega = 2.*pi/lambda*c0;    % the angular frequency

k = 2*pi/lambda*nm;
M = numel(an);
n = (1:M)';

E0 = 1;
En = 1j.^n.*E0.*(2*n+1)./(n.*(n+1));

%% Get internal field expansion coefficients
if strtfd
    [fn, gn, vn, wn] = expcoeff_mie_strat_int(an, bn, x, m);
else %if strtfd
    [fn, gn] = expcoeff_mie_int(x, m, numel(an));
    vn = zeros(size(fn));
    wn = zeros(size(fn));
end %if strtfd

E = zeros([numel(xf), 3]);
H = zeros([numel(xf), 3]);

for ic=1:numel(xf)
    
    % Get spherical coordinates
    [theta, phi, rad] = xcart2sph(xf(ic), yf(ic), zf(ic));
    
    % Calculate angle dependent functions
    [pin, taun] = angdepfun_mie(theta, n);
    
    % calculate coordinate transformation matrix
    R = getTransformationMatrix(phi, theta);
    
    % Initialize arrays
    M_o1n = zeros(M,3);
    M_e1n = zeros(M,3);
    N_o1n = zeros(M,3);
    N_e1n = zeros(M,3);
    
    Edum = zeros(M,3);
    Hdum = zeros(M,3);
    
    if rad > r(end)
        % external fields
        
        % Auxiliary parameters
        rho = k*rad;
        h = sbesselh(n, rho);
        dxi = dricbesh(n, rho);
        
        M_o1n(:,2) = cos(phi)*pin.*h;
        M_o1n(:,3) = -sin(phi)*taun.*h;
        M_e1n(:,2) = -sin(phi)*pin.*h;
        M_e1n(:,3) = -cos(phi)*taun.*h;
        N_o1n(:,1) = sin(phi)*n.*(n + 1)*sin(theta).*pin.*h/rho;
        N_o1n(:,2) = sin(phi).*taun.*dxi/rho;
        N_o1n(:,3) = cos(phi).*pin.*dxi/rho;
        N_e1n(:,1) = cos(phi)*n.*(n + 1)*sin(theta).*pin.*h/rho;
        N_e1n(:,2) = cos(phi).*taun.*dxi/rho;
        N_e1n(:,3) = -sin(phi).*pin.*dxi/rho;
        
        for d=1:3
            Edum(:,d) = En.*(1.j*an.'.*N_e1n(:,d) - bn.'.*M_o1n(:,d));
            Hdum(:,d) = En.*(1.j*bn.'.*N_o1n(:,d) + an.'.*M_e1n(:,d));
        end %for d=1:3
        
        E(ic,:) = R*sum(Edum).';
        H(ic,:) = k/omega/mue0*R*sum(Hdum).';
        
        if tf_flag
            % total field
            E(ic,1) = E(ic,1) + exp(1.j*k*zf(ic));
            H(ic,2) = H(ic,2) + nm/c0/mue0*exp(1.j*k*zf(ic));
        end %if tf_flag
        
        if ~cc_flag
            E(ic,:) = R'*E(ic,:).';
            H(ic,:) = R'*H(ic,:).';
        end %if ~cc_flag
        
    else %if rad > r
        % internal fields
        ilay = max([find(r >=rad, 1, 'first') , find(r < rad, 1, 'last')]);
        
        % Auxiliary parameters
        rho = k1(ilay)*rad;
        h = sbesselj(n, rho);
        dxi = dricbesj(n, rho);
        
        M_o1n(:,2) = cos(phi)*pin.*h;
        M_o1n(:,3) = -sin(phi)*taun.*h;
        M_e1n(:,2) = -sin(phi)*pin.*h;
        M_e1n(:,3) = -cos(phi)*taun.*h;
        N_o1n(:,1) = sin(phi)*n.*(n + 1)*sin(theta).*pin.*h/rho;
        N_o1n(:,2) = sin(phi).*taun.*dxi/rho;
        N_o1n(:,3) = cos(phi).*pin.*dxi/rho;
        N_e1n(:,1) = cos(phi)*n.*(n + 1)*sin(theta).*pin.*h/rho;
        N_e1n(:,2) = cos(phi).*taun.*dxi/rho;
        N_e1n(:,3) = -sin(phi).*pin.*dxi/rho;
        
        for d=1:3
            Edum(:,d) = En.*(fn(:,ilay).*M_o1n(:,d) - ...
                1.j*gn(:,ilay).*N_e1n(:,d));
            Hdum(:,d) = En.*(gn(:,ilay).*M_e1n(:,d) + ...
                1.j*fn(:,ilay).*N_o1n(:,d));
        end %for d=1:3
        
        h = sbessely(n, rho);
        dxi = -dricbesy(n, rho);
        
        M_o1n(:,2) = cos(phi)*pin.*h;
        M_o1n(:,3) = -sin(phi)*taun.*h;
        M_e1n(:,2) = -sin(phi)*pin.*h;
        M_e1n(:,3) = -cos(phi)*taun.*h;
        N_o1n(:,1) = sin(phi)*n.*(n + 1)*sin(theta).*pin.*h/rho;
        N_o1n(:,2) = sin(phi).*taun.*dxi/rho;
        N_o1n(:,3) = cos(phi).*pin.*dxi/rho;
        N_e1n(:,1) = cos(phi)*n.*(n + 1)*sin(theta).*pin.*h/rho;
        N_e1n(:,2) = cos(phi).*taun.*dxi/rho;
        N_e1n(:,3) = -sin(phi).*pin.*dxi/rho;
        
        for d=1:3
            Edum(:,d) = Edum(:,d) + En.*(vn(:,ilay).*M_o1n(:,d) - ...
                1.j*wn(:,ilay).*N_e1n(:,d));
            Hdum(:,d) = Hdum(:,d) + En.*(wn(:,ilay).*M_e1n(:,d) + ...
                1.j*vn(:,ilay).*N_o1n(:,d));
        end %for d=1:3
        
        E(ic,:) = R*sum(Edum).';
        H(ic,:) = -k1(ilay)/omega/mue0*R*sum(Hdum).';
        
        if ~tf_flag
            % total field
            E(ic,1) = E(ic,1) - exp(1.j*k*zf(ic));
            H(ic,2) = H(ic,2) - nm/c0/mue0*exp(1.j*k*zf(ic));
        end %if ~tf_flag
        
        if ~cc_flag
            E(ic,:) = R'*E(ic,:).';
            H(ic,:) = R'*H(ic,:).';
        end %if ~cc_flag
    end %if rad > r
    
end %for ic=1:numel(xc)

E = squeeze(reshape(E,[size(xf), 3]));
H = squeeze(reshape(H,[size(xf), 3]));

end

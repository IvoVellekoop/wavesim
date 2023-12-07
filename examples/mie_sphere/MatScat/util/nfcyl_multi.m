function [ Ep, En, Hp, Hn ] = nfcyl_multi( ann, bnp, Ann, Bnp, r, ...
    ns, nm, lambda, xc, yc, xf, yf, iphi, tf_flag, cc_flag )
%NFCYL_MULTI Calculates the near field solution for given expansion
%   coefficients of many cylindrical particles.
%
%   [EP,EN,HP,HN] = NFCYL_MULTI(ann,bnp,Ann,Bnp,R,NS,NM,LAMBDA,XC,YC,XF,YF,
%   IPHI,TF_FLAG,CC_FLAG) calculates the near field E and H for the 
%   scattering of a monochromatic electromagnetic wave at perpendicular 
%   incidence by many cylindrical particles with given scattering expansion
%   coefficients ann and bnp and internal expansion coefficients Ann and
%   Bnp. The solution is for incoming electromagnetic waves polarized 
%   perpendicular (EP,HP) and normal (EN, HN) according to the cylinder 
%   axes. The cylinders have radius R and (complex) refractive index NS and
%   are located at coordinates given by XC and YC. The incident wave has a 
%   wavelength of LAMBDA and is propagating at an incident angle IPHI (in
%   radians) .The near field solution is calculated for 3D cartesian 
%   coordinates given in arrays XF and YF.
%
%   If TF_FLAG is set total fields are returned, otherwise scattered fields
%   are calculated. IF CC_FLAG is set resulting vectors are in cartesian
%   coordinate vector basis, in cylindrical coordinate basis otherwise.
%
%   The calculation is based on the derivation of Lee [2]:
%   [2] Lee, S.-C., Dependent scattering of an obliquely incident plane
%       wave by a collection of parallel cylinders. J. Appl. Phys. 68(10),
%       1990.
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Initialize parameters
if ~exist('iphi', 'var')
    iphi = 0.;
end %if ~exist('iphi', 'var')

if ~exist('tf_flag', 'var')
    tf_flag = false;
end %if ~exist('tf_flag', 'var')

if ~exist('cc_flag', 'var')
    cc_flag = true;
end %if ~exist('cc_flag', 'var')

MatScat_const;

hk = 2;                 % allways use Hankel function of the 2nd kind

k = 2*pi/lambda*nm;     % wavenumber in medium nm
k1 = 2*pi/lambda*ns;    % wavenumber in medium ns
N = size(ann, 1);       % total number of cylidners

%% Get truncation number
MM = size(ann,2);
M = (MM - 1)/2;
n = (-M:M)';
ng = meshgrid(n, 1:N)';

%% Define field arrays
Ep = zeros([numel(xf), 3]);
En = zeros([numel(xf), 3]);
Hp = zeros([numel(xf), 3]);
Hn = zeros([numel(xf), 3]);

for ic=1:numel(xf)
    
    % Calculate distance and angle matrices
    dist = zeros(1,N);
    gamkp = zeros(1,N);
    for i=1:N
        rf = [xf(ic), yf(ic)];
        rc = [xc(i), yc(i)];
        rv = rf - rc;
        dist(i) = norm(rv);
        gamkp(i) = acos(dot(rv, [1, 0])/dist(i));
        if isnan(gamkp(i))
            gamkp(i) = 0;
        end %if isnan(gamkp(i))
        if rv(2) < 0
            gamkp(i) = -gamkp(i);
        end %if rv(2) < 0
    end %for i=1:N
    
    % get polar coordinates
    [phi, ~, ~] = cart2pol(xf(ic), yf(ic), 0);
    
    % calculate coordinate transformation matrix
    R = getTransformationMatrix(phi);
    
    if all(dist >= r)
        % scattered field
        
        % auxiliary variables
        expdum = (-1.j).^ng.*exp(1.j*n*gamkp);
        [nu,z] = meshgrid(n', k*dist');
        bh = besselh(nu, hk, z).';
        dbh = dbesselh(nu, hk, z).';
        rhoj = expdum.*(1.j.*n*(1./dist)).*bh;
        phij = k*expdum.*dbh;
        
        Erhoj = -rhoj.*ann.';
        Ephij = phij.*ann.';
        
        Hrhoj = rhoj.*bnp.';
        Hphij = -phij.*bnp.';
        
        Erho = sum(sum(Erhoj).*cos(gamkp-phi) - ...
            sum(Ephij).*sin(gamkp-phi));
        Ephi = sum(sum(Erhoj).*sin(gamkp-phi) + ...
            sum(Ephij).*cos(gamkp-phi));
        Ez = k*1.j*sum(sum(-expdum.*bh.*bnp.'));
        
        Hrho = nm*sum(sum(Hrhoj).*cos(gamkp-phi) - ...
            sum(Hphij).*sin(gamkp-phi));
        Hphi = nm*sum(sum(Hrhoj).*sin(gamkp-phi) + ...
            sum(Hphij).*cos(gamkp-phi));
        Hz = k*1.j*sum(sum(-expdum.*nm.*bh.*ann.'));
        
        Ep(ic,:) = R*[0, 0, Ez].'/1.j./k;
        En(ic,:) = R*[Erho, Ephi, 0].'/1.j./k;
        Hp(ic,:) = R*[Hrho, Hphi, 0].'/1.j./k/c0/mue0;
        Hn(ic,:) = R*[0., 0., Hz].'/1.j./k/c0/mue0;
               
        if tf_flag
            % total field
            epsc = exp(-1.i*k*(xf(ic)*cos(iphi) - yf(ic)*sin(iphi)));            
            Ep(ic,3) = Ep(ic,3) + epsc;
            Hp(ic,1) = Hp(ic,1) - nm/c0/mue0*epsc*sin(iphi);
            Hp(ic,2) = Hp(ic,2) - nm/c0/mue0*epsc*cos(iphi);
            En(ic,1) = En(ic,1) + epsc*sin(iphi);
            En(ic,2) = En(ic,2) + epsc*cos(iphi);
            Hn(ic,3) = Hn(ic,3) + nm/c0/mue0*epsc;
        end %if tf_flag
        
        if ~cc_flag
            Ep(ic,:) = R'*Ep(ic,:).';
            En(ic,:) = R'*En(ic,:).';
            Hp(ic,:) = R'*Hp(ic,:).';
            Hn(ic,:) = R'*Hn(ic,:).';
        end %if ~cc_flag
        
    else %if all(dist >= r)
        % internal total field
        i = find(dist < r);
        
        % auxiliary variables                
        expdum = (-1.j).^n.*exp(1.j*n*gamkp(i));
        bj = besselj(n, k1*dist(i));
        rhoj = expdum.*(1.j.*n/dist(i)).*bj;
        phij = k1*expdum.*dbesselj(n, k1*dist(i));
        
        Erhoj = rhoj.*Ann(i,:).';
        Ephij = -phij.*Ann(i,:).';
        
        Hrhoj = -rhoj.*Bnp(i,:).';
        Hphij = phij.*Bnp(i,:).';
        
        Erho = sum(Erhoj).*cos(gamkp(i)-phi) - ...
            sum(Ephij).*sin(gamkp(i)-phi);
        Ephi = sum(Erhoj).*sin(gamkp(i)-phi) + ...
            sum(Ephij).*cos(gamkp(i)-phi);
        Ez = sum(expdum.*Bnp(i,:).'.*k1.*1.j.*bj);
        
        Hrho = ns*(sum(Hrhoj).*cos(gamkp(i)-phi) - ...
            sum(Hphij).*sin(gamkp(i)-phi));
        Hphi = ns*(sum(Hrhoj).*sin(gamkp(i)-phi) + ...
            sum(Hphij).*cos(gamkp(i)-phi));
        Hz = sum(expdum.*Ann(i,:).'.*k1.*1.j.*ns.*bj);
        
        Ep(ic,:) = R*[0, 0, Ez].'/1.j./k;
        Hp(ic,:) = R*[Hrho, Hphi, 0].'/1.j./k/c0/mue0;
        En(ic,:) = R*[Erho, Ephi, 0.].'/1.j./k;
        Hn(ic,:) = R*[0., 0., Hz].'/1.j./k/c0/mue0;
        
        if ~tf_flag
            % internal scattered field
            epsc = exp(-1.i*k*(xf(ic)*cos(-iphi) + yf(ic)*sin(-iphi)));            
            Ep(ic,3) = Ep(ic,3) - epsc;
            Hp(ic,1) = Hp(ic,1) + nm/c0/mue0*epsc*sin(iphi);
            Hp(ic,2) = Hp(ic,2) + nm/c0/mue0*epsc*cos(iphi);
            En(ic,1) = En(ic,1) - epsc*sin(iphi);
            En(ic,2) = En(ic,2) - epsc*cos(iphi);
            Hn(ic,3) = Hn(ic,3) - nm/c0/mue0*epsc;
        end %if ~tf_flag
        
        if ~cc_flag
            Ep(ic,:) = R'*Ep(ic,:).';
            En(ic,:) = R'*En(ic,:).';
            Hp(ic,:) = R'*Hp(ic,:).';
            Hn(ic,:) = R'*Hn(ic,:).';
        end %if ~cc_flag
    end %if all(dist >= r)
end %for icoord=1:length(coords)

Ep = reshape(Ep,[size(xf), 3]);
En = reshape(En,[size(xf), 3]);
Hp = reshape(Hp,[size(xf), 3]);
Hn = reshape(Hn,[size(xf), 3]);

end

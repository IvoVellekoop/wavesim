function [ Sp, Sn ] = nfcyl_multi_Poynting( ann, bnp, Ann, Bnp, r, ns, ...
nm, lambda, xc, yc, xf, yf, tf_flag, cc_flag )
%NFCYL_MULTI_POYNTING Calculates the Poynting vector for given expansion
%   coefficients of many cylindrical particles.
%
%   [SP,SN] = NFCYL_MULTI_POYNTING(ann,bnp,Ann,Bnp,R,NS,NM,LAMBDA,XC,YC,XF,
%   YF,IPHI,TF_FLAG,CC_FLAG) calculates the Poynting Vector for the 
%   scattering of a monochromatic electromagnetic wave at perpendicular 
%   incidence by many cylindrical particles with given scattering expansion
%   coefficients ann and bnp and internal expansion coefficients Ann and
%   Bnp. The solution is for incoming electromagnetic waves polarized 
%   perpendicular SP and normal SN according to the cylinder axes. The 
%   cylinders have radius R and (complex) refractive index NS and are 
%   located at coordinates given by XC and YC. The incident wave has a 
%   wavelength of LAMBDA. The near field solution is calculated for 3D 
%   cartesian coordinates given in arrays XF and YF.
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
if ~exist('tf_flag', 'var')
    tf_flag = false;
end %if ~exist('tf_flag', 'var')

if ~exist('cc_flag', 'var')
    cc_flag = true;
end %if ~exist('cc_flag', 'var')

Mielab_globals;

hk = 2;                 % allways use Hankel function of the 2nd kind

k = 2*pi/lambda*nm;     % wavenumber in medium nm
k1 = 2*pi/lambda*ns;    % wavenumber in medium ns
N = size(ann, 1);       % total number of cylidners

%% Get truncation number
MM = size(ann,2);
M = (MM - 1)/2;
n = (-M:M)';
ns = n*ones(size(n))';

%% Define field arrays
Sp = zeros([numel(xf), 3]);
Sn = zeros([numel(xf), 3]);

bnp = bnp.';
ann = ann.';
Bnp = Bnp.';
Ann = Ann.';

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
    if cc_flag
        R = getRotationMatrix(phi);
    else %if cc_flag
        R = getRotationMatrix(0);
    end %if cc_flag
    
    if all(dist > r)
        if ~tf_flag
            Sprho = 0.;
            Spphi = 0.;
            Snrho = 0.;
            Snphi = 0.;
            for i=1:N
                phii = phi-gamkp(i);
                for j=1:N
                    phij = phi-gamkp(j);
                    expdum = (-1.j).^(ns-ns').*...
                        exp(1.j*(ns*gamkp(i) - ns'*gamkp(j)));
                    [nu,z] = meshgrid(ns, k*dist(i));
                    bhi = besselh(nu, hk, z);
                    bhj = besselh(nu', hk, z);
                    dbhi = dbesselh(nu, hk, z);
                    dbhj = dbesselh(nu', hk, z);
                    bn2 = bnp(:,i)*bnp(:,j)';
                    an2 = ann(:,i)*ann(:,j)';
                    Sprho = Sprho + sum(sum(expdum.*...
                        (bhi.*conj(dbhj).*cos(phij) - ...
                        bhi.*conj(bhj).*1.j.*ns'/k/dist(j).*sin(phij)).*bn2));
                    Spphi = Spphi + sum(sum(expdum.*...
                        (bhi.*conj(bhj).*ns'/k/dist(j).*cos(phij) - ...
                        bhi.*conj(dbhj).*1.j.*sin(phij)).*bn2));
                    Snrho = Snrho + sum(sum(expdum.*...
                        (dbhi.*conj(bhj).*cos(phii) + ...
                        bhi.*conj(bhj).*1.j.*ns/k/dist(i).*sin(phii)).*an2));
                    Snphi = Snphi + sum(sum(expdum.*...
                        (bhi.*conj(bhj).*ns/k/dist(i).*cos(phii) + ...
                        dbhi.*conj(bhj).*1.j.*sin(phii)).*an2));
                end %for j=1:N
            end %for i=1:N
            
            Sprho = -real(1.j*Sprho);
            Spphi = -real(Spphi);
            Snrho = real(1.j*Snrho);
            Snphi = -real(Snphi);
            
            Sp(ic,:) = R*[Sprho, Spphi, 0.].';
            Sn(ic,:) = R*[Snrho, Snphi, 0.].';
        end %if ~tf_flag
    else %if all(dist >= r)
        if tf_flag
            j = find(dist < r);
            
            phij = phi-gamkp(j);
            expdum = (-1.j).^(ns-ns').*...
                exp(1.j*((ns-ns').*gamkp(j)));
            bhn = besselj(ns, k1*dist(j));
            bhs = besselj(ns', k1*dist(j));
            dbhn = dbesselj(ns, k1*dist(j));
            dbhs = dbesselj(ns', k1*dist(j));
            bn2 = Bnp(:,j)*Bnp(:,j)';
            an2 = Ann(:,j)*Ann(:,j)';
            Sprho = sum(sum(expdum.*...
                (bhn.*conj(dbhs).*cos(phij) + ...
                bhn.*conj(bhs).*1.j.*ns'/k1/dist(j).*sin(phij)).*bn2));
            Spphi = sum(sum(expdum.*...
                (bhn.*conj(bhs).*ns'/k1/dist(j).*cos(phij) + ...
                bhn.*conj(dbhs).*1.j.*sin(phij)).*bn2));
            Snrho = sum(sum(expdum.*...
                (dbhn.*conj(bhs).*cos(phij) + ...
                bhn.*conj(bhs).*1.j.*ns/k1/dist(j).*sin(phij)).*an2));
            Snphi = sum(sum(expdum.*...
                (bhn.*conj(bhs).*ns/k1/dist(j).*cos(phij) + ...
                dbhn.*conj(bhs).*1.j.*sin(phij)).*an2));
            
            Sprho = -real(1.j*Sprho);
            Spphi = -real(Spphi);
            Snrho = real(1.j*Snrho);
            Snphi = -real(Snphi);
            
            Sp(ic,:) = R*[Sprho, Spphi, 0.].';
            Sn(ic,:) = R*[Snrho, Snphi, 0.].';
        end %if tf_flag
    end %if all(dist >= r)
end %for icoord=1:length(coords)


Sp = reshape(Sp,[size(xf), 3]);
Sn = reshape(Sn,[size(xf), 3]);

end

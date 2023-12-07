function [fnp, gnn, vnp, wnn] = expcoeff_cyl_strat_int( ann, bnp, r, k, m )
%EXPCOEFF_CYL_STRAT_INT Calculates the internal expansion coefficients for 
%   the scattering by a stratified infinite cylinder at perpendicular
%   incidence.
%
%   [FNP,GNN,VNP,WNN] = EXPCOEFF_CYL_STRAT_INT(ANN,BNP,R,K,M) calculates 
%   the expansion coefficients (FNP,GNN,VNP,WNN) for a stratified cylinder 
%   with size parameters X and relative refractive indices M at 
%   perpendicular incidence with wavenumber K. 
%
%   The calculation is based on the book of Kerker [5]:
%   [5] Kerker, M., The scattering of light and other electromagnetic
%       radiation, Academic Press, 1969
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)


N = numel(r);   % Number of layers
MM = numel(ann);
M = (MM - 1)/2;

%% Calculate expansion coefficients
fnp = zeros(MM,N);
gnn = zeros(MM,N);
vnp = zeros(MM,N);
wnn = zeros(MM,N);

for ilay=N:-1:1
    if ilay == N
        mfrac = 1/m(ilay);
        rho_lp1 = k*r(ilay);
        lp1 = k;
    else
        mfrac =  m(ilay+1)/m(ilay); 
        lp1 = k*m(ilay+1);
        rho_lp1 = lp1*r(ilay);       
    end
    rho_l = k*m(ilay)*r(ilay);
    for n=-M:M
        if ilay == N
            u_lp1 = besselj(n, rho_lp1) - bnp(n+M+1)*besselh(n, rho_lp1);
            u_lp1_dr = lp1*(dbesselj(n, rho_lp1) - ...
                bnp(n+M+1)*dbesselh(n, rho_lp1));
            v_lp1 = besselj(n, rho_lp1) - ann(n+M+1)*besselh(n, rho_lp1);
            v_lp1_dr = lp1*(dbesselj(n, rho_lp1) - ...
                ann(n+M+1)*dbesselh(n, rho_lp1));
        else
            u_lp1 = fnp(n+M+1,ilay+1)*besselj(n, rho_lp1) - ...
                vnp(n+M+1,ilay+1)*besselh(n, rho_lp1);
            u_lp1_dr = lp1*(fnp(n+M+1,ilay+1)*dbesselj(n, rho_lp1) - ...
                vnp(n+M+1,ilay+1)*dbesselh(n, rho_lp1));
            v_lp1 = gnn(n+M+1,ilay+1)*besselj(n, rho_lp1) - ...
                wnn(n+M+1,ilay+1)*besselh(n, rho_lp1);
            v_lp1_dr = lp1*(gnn(n+M+1,ilay+1)*dbesselj(n, rho_lp1) - ...
                wnn(n+M+1,ilay+1)*dbesselh(n, rho_lp1));
        end
        if ilay > 1
            denom = dbesselh(n, rho_l)*besselj(n, rho_l) - ...
                besselh(n, rho_l)*dbesselj(n, rho_l);
            
            vnp(n+M+1,ilay) = mfrac/m(ilay)/k*(k*m(ilay)*...
                dbesselj(n, rho_l)*u_lp1 - ...
                besselj(n, rho_l)*u_lp1_dr)/denom;
            fnp(n+M+1,ilay) = mfrac/m(ilay)/k*(k*m(ilay)*...
                dbesselh(n, rho_l)*u_lp1 - ...
                besselh(n, rho_l)*u_lp1_dr)/denom;
            wnn(n+M+1,ilay) = 1/m(ilay)/k*(lp1*mfrac*...
                dbesselj(n, rho_l)*v_lp1 - ...
                besselj(n, rho_l)*v_lp1_dr)/denom;
            gnn(n+M+1,ilay) = 1/m(ilay)/k*(lp1*mfrac*...
                dbesselh(n, rho_l)*v_lp1 - ...
                besselh(n, rho_l)*v_lp1_dr)/denom;
        else
            fnp(n+M+1,ilay) = mfrac*u_lp1/besselj(n, rho_l);
            gnn(n+M+1,ilay) = mfrac.^2*v_lp1/besselj(n, rho_l);
        end %ilay
    end %for n=1:M
end %for ilay=N:-1:1

end


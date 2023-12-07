function [fn, gn, vn, wn] = expcoeff_mie_strat_int( an, bn, x, m )
%EXPCOEFF_MIE_STRAT_INT Calculates the expansion coefficients of the 
%   internal fields for a stratified sphere.
%
%   [FN,GN,VN,WN] = EXPCOEFF_MIE_STRAT(AN,BN,X,M) calculates the internal 
%   fields expansion coefficients FN, GN, VN and WN for a stratified sphere
%   with size parameters X and relative refractive indices M. The
%   corresponding scattered field expansion coefficients are given in AN
%   and BN.
%
%   The calculation is based on the book of Bohren and Huffman [1]:
%   [1] Bohren, C. F. and Huffman, D. R., Absorption and scattering of
%       light by small particles, Wiley-Interscience, New York, 1998.
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Calculate truncation number
M = numel(an);
N = numel(x);

%% Calculate expansion coefficients
fn = zeros(M,N);
gn = zeros(M,N);
vn = zeros(M,N);
wn = zeros(M,N);

for ilay=N:-1:1
    for n=1:M
        if ilay == N
            Sx = ricbesj(n, x(ilay));
            dSx = dricbesj(n, x(ilay));
            Smx = ricbesj(n, m(ilay)*x(ilay));
            dSmx = dricbesj(n, m(ilay)*x(ilay));
            xix = ricbesh(n, x(ilay));
            dxix = dricbesh(n, x(ilay));
            cmx = ricbesy(n, m(ilay)*x(ilay));
            dcmx = dricbesy(n, m(ilay)*x(ilay));
            denom = cmx*dSmx - dcmx*Smx;
            fn(n,ilay) = ((Sx - bn(n)*xix)*dcmx*m(ilay) + ...
                (bn(n)*dxix - dSx)*cmx)/-denom;
            gn(n,ilay) = ((dSx - an(n)*dxix)*cmx*m(ilay) + ...
                (an(n)*xix - Sx)*dcmx)/denom;
            vn(n,ilay) = ((dSx - bn(n)*dxix)*Smx + ...
                (bn(n)*xix - Sx)*dSmx*m(ilay))/denom;
            wn(n,ilay) = ((Sx - an(n)*xix)*dSmx + ...
                (an(n)*dxix - dSx)*Smx*m(ilay))/-denom;
        elseif ilay == 1
            Smp1x = ricbesj(n, m(ilay+1)*x(ilay));
            Smx = ricbesj(n, m(ilay)*x(ilay));
            cmp1x = ricbesy(n, m(ilay+1)*x(ilay));
            fn(n,ilay) = (fn(n,ilay+1)*m(ilay)*Smp1x - ...
                vn(n,ilay+1)*m(ilay)*cmp1x)/m(ilay+1)/Smx;
            gn(n,ilay) = (gn(n,ilay+1)*Smp1x - ...
                wn(n,ilay+1)*cmp1x)/Smx;
            
            % dSmx = dricbesj(n, m(ilay)*x(ilay));
            % dSmp1x = dricbesj(n, m(ilay+1)*x(ilay));
            % dcmp1x = dricbesy(n, m(ilay+1)*x(ilay));
            % fn(n,ilay) = (fn(n,ilay+1)*dSmp1x - ...
            %     vn(n,ilay+1)*dcmp1x)/dSmx;
            % gn(n,ilay) = (gn(n,ilay+1)*m(ilay)*dSmp1x - ...
            %     wn(n,ilay+1)*m(ilay)*dcmp1x)/m(ilay+1)/dSmx;
            
        else %ilay
            Smp1x = ricbesj(n, m(ilay+1)*x(ilay));
            dSmp1x = dricbesj(n, m(ilay+1)*x(ilay));
            Smx = ricbesj(n, m(ilay)*x(ilay));
            dSmx = dricbesj(n, m(ilay)*x(ilay));
            cmp1x = ricbesy(n, m(ilay+1)*x(ilay));
            dcmp1x = dricbesy(n, m(ilay+1)*x(ilay));
            cmx = ricbesy(n, m(ilay)*x(ilay));
            dcmx = dricbesy(n, m(ilay)*x(ilay));
            denom = m(ilay+1)*(Smx*dcmx - dSmx*cmx);
            fn(n,ilay) = (fn(n,ilay+1)*(m(ilay)*Smp1x*dcmx - ...
                m(ilay+1)*dSmp1x*cmx) + ...
                vn(n,ilay+1)*(m(ilay+1)*dcmp1x*cmx - ...
                m(ilay)*cmp1x*dcmx))/denom;
            gn(n,ilay) = (gn(n,ilay+1)*(m(ilay+1)*Smp1x*dcmx - ...
                m(ilay)*dSmp1x*cmx) + ...
                wn(n,ilay+1)*(m(ilay)*dcmp1x*cmx - ...
                m(ilay+1)*cmp1x*dcmx))/denom;
            vn(n,ilay) = (fn(n,ilay+1)*(m(ilay)*Smp1x*dSmx - ...
                m(ilay+1)*dSmp1x*Smx) + ...
                vn(n,ilay+1)*(m(ilay+1)*dcmp1x*Smx - ...
                m(ilay)*cmp1x*dSmx))/denom;
            wn(n,ilay) = (gn(n,ilay+1)*(m(ilay+1)*Smp1x*dSmx - ...
                m(ilay)*dSmp1x*Smx) + ...
                wn(n,ilay+1)*(m(ilay)*dcmp1x*Smx - ...
                m(ilay+1)*cmp1x*dSmx))/denom;
        end %ilay
    end %for n=1:M
end %for ilay=N:-1:1

end

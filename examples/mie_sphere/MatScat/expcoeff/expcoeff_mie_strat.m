function [an, bn] = expcoeff_mie_strat( x, m, conv )
%EXPCOEFF_MIE_STRAT Calculates the expansion coefficients for the
%   scattering by a stratified sphere.
%
%   [AN,BN] = EXPCOEFF_MIE_STRAT(X,M,CONV) calculates the expansion 
%   coefficients AN and BN for a stratified sphere with size parameters and
%   relative refractive indices for the layers given in lists X and M.
%   The convergence criteria is multiplied by a factor of CONV.
%
%   The calculation is based on the book of Kerker [5]:
%   [5] Kerker M., The scattering of light and other electromagnetic
%       radiation, Academic Press, New York, 1969.
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

%% Initialize parameters
if ~exist('conv', 'var')
    conv = 1;
end %if ~exist('conv', 'var')

%% Calculate truncation number
M = ceil(conv*(x(end) + 4*(x(end)^(1/3)) + 2));
N = numel(x);

%% Initialize arrays
an = zeros(1,M);
bn = zeros(1,M);

if N == 2
    
    n=1:M;
    
    Sm1x = ricbesj(n, m(1)*x(1));
    dSm1x = dricbesj(n, m(1)*x(1));
    Sm2x = ricbesj(n, m(2)*x(1));
    dSm2x = dricbesj(n, m(2)*x(1));
    cm2x = ricbesy(n, m(2)*x(1));
    dcm2x = dricbesy(n, m(2)*x(1));
    
    An = (m(2)*Sm2x.*dSm1x - m(1)*dSm2x.*Sm1x)./...
        (m(2)*cm2x.*dSm1x - m(1)*dcm2x.*Sm1x);
    Bn = (m(2)*Sm1x.*dSm2x - m(1)*dSm1x.*Sm2x)./...
        (m(2)*dcm2x.*Sm1x - m(1)*cm2x.*dSm1x);
    
    Sy = ricbesj(n, x(2));
    dSy = dricbesj(n, x(2));
    Sm2y = ricbesj(n, m(2)*x(2));
    dSm2y = dricbesj(n, m(2)*x(2));
    cm2y = ricbesy(n, m(2)*x(2));
    dcm2y = dricbesy(n, m(2)*x(2));
    xiy = ricbesh(n, x(2));
    dxiy = dricbesh(n, x(2));
    an = (Sy.*(dSm2y - An.*dcm2y) - m(2)*dSy.*(Sm2y - An.*cm2y))./...
        (xiy.*(dSm2y - An.*dcm2y) - m(2)*dxiy.*(Sm2y - An.*cm2y));
    bn = (m(2)*Sy.*(dSm2y - Bn.*dcm2y) - dSy.*(Sm2y - Bn.*cm2y))./...
        (m(2)*xiy.*(dSm2y - Bn.*dcm2y) - dxiy.*(Sm2y - Bn.*cm2y));
    
else %if N= 2
    
    A = zeros(2*N,2*N);
    C = zeros(2*N,2*N);
    
    for n=1:M
        % Setup the equation system
        for i=1:2*N
            for j=1:2*N
                p = floor(j/2) + 1;
                q = floor((i + 1)/2);
                if (p-q == 0) || (p-q == 1)
                    if mod(i,2) == 1
                        if (j < 2*N) && ((j == 1) || (mod(j,2) == 0))
                            A(i,j) = dricbesj(n, m(p)*x(q));
                        else if (mod(j,2) == 1)
                                A(i,j) = dricbesy(n, m(p)*x(q));
                            else
                                A(i,j) = dricbesj(n, x(q));
                            end %if (mod(i,2) == 0)
                        end %if (j < 2*N) && ((j == 1) || (mod(i,2) == 0))
                        if j ~= 2*N
                            C(i,j) = m(p)*A(i,j);
                        else %if j ~= 2*N
                            C(i,j) = A(i,j);
                        end %if j ~= 2*N
                    else %if mod(i,2) == 1
                        if (j < 2*N) && ((j == 1) || (mod(j,2) == 0))
                            C(i,j) = ricbesj(n, m(p)*x(q));
                        else if (mod(j,2) == 1)
                                C(i,j) = ricbesy(n, m(p)*x(q));
                            else
                                C(i,j) = ricbesj(n, x(q));
                            end %if (mod(i,2) == 0)
                        end %if (j < 2*N) && ((j == 1) || (mod(i,2) == 0))
                        if j ~= 2*N
                            A(i,j) = m(p)*C(i,j);
                        else %if j ~= 2*N
                            A(i,j) = C(i,j);
                        end %if j ~= 2*N
                    end %if mod(i,2) == 1
                end %if (p-q == 0) || (p-q == 1)
            end %for j=1:2*N
        end %for i=1:2*N
        
        % Solve the equation system
        B = A;
        B(end-1,end) = dricbesh(n, x(end));
        B(end,end) = ricbesh(n, x(end));
        an(n) = det(A)/det(B);
        
        D = C;
        D(end-1,end) = dricbesh(n, x(end));
        D(end,end) = ricbesh(n, x(end));
        bn(n) = det(C)/det(D);
        
    end %for n=-M:M
end %if N= 2
end


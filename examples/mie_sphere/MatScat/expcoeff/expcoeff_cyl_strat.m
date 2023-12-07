function [ann, bnp] = expcoeff_cyl_strat( x, m, hk, conv )
%EXPCOEFF_CYL_STRAT Calculates the expansion coefficients for the 
%   scattering by a stratified infinite cylinder at perpendicular 
%   incidence.
%
%   [ANN,BNP] = EXPCOEFF_CYL_STRAT(X,M,HK,CONV) calculates the expansion 
%   coefficients (ANN, BNP) for a stratified cylinder with size parameters 
%   X and relative refractive indices M. The Hankel function H_N^(HK) is 
%   used in the computation. The convergence criteria is multiplied by a 
%   factor of CONV.
%
%   The calculation is based on the book of Kerker [5]:
%   [5] Kerker, M., The scattering of light and other electromagnetic 
%       radiation, Academic Press, 1969 
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)


%% Initialize parameters
if ~exist('hk', 'var')
    hk = 1;
end %if ~exist('hk', 'var')

if ~exist('conv', 'var')
    conv = 1;
end %if ~exist('conv', 'var')

N = numel(x);   % Number of layers

%% Calculate truncation number
M = ceil(conv*ceil(x(end) + 4*(x(end)^(1/3)) + 2));

%% Calculate expansion coefficients
ann = zeros(1,M+1);
bnp = zeros(1,M+1);

A = zeros(2*N,2*N);
C = zeros(2*N,2*N);
for n=0:M
    for i=1:2*N
        for j=1:2*N
            p = floor(j/2) + 1;
            q = floor((i + 1)/2);
            if (p-q == 0) || (p-q == 1)
                if mod(i,2) == 1
                    if (j < 2*N) && ((j == 1) || (mod(j,2) == 0))
                        A(i,j) = dbesselj(n, m(p)*x(q));
                    else if (mod(j,2) == 1)
                            A(i,j) = dbessely(n, m(p)*x(q));
                        else
                            A(i,j) = dbesselj(n, x(q));
                        end %if (mod(i,2) == 0)
                    end %if (j < 2*N) && ((j == 1) || (mod(i,2) == 0))
                    if j ~= 2*N                       
                        C(i,j) = m(p)*A(i,j);
                    else %if j ~= 2*N                       
                        C(i,j) = A(i,j);
                    end %if j ~= 2*N
                else %if mod(i,2) == 1
                    if (j < 2*N) && ((j == 1) || (mod(j,2) == 0))
                        C(i,j) = besselj(n, m(p)*x(q));
                    else if (mod(j,2) == 1)
                            C(i,j) = bessely(n, m(p)*x(q));
                        else
                            C(i,j) = besselj(n, x(q));
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
  
    B = A;
    B(end-1,end) = dbesselh(n, hk, x(end));
    B(end,end) = besselh(n, hk, x(end));
    ann(n+1) = det(A)/det(B);
    
    D = C;
    D(end-1,end) = dbesselh(n, hk, x(end));
    D(end,end) = besselh(n, hk, x(end));
    bnp(n+1) = det(C)/det(D);

end %for n=-M:M

end


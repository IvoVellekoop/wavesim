function [pin, taun] = angdepfun_mie( theta, n )
%ANGDEPFUN_MIE Calculates the angle dependent functions for the Mie theory.
%
%   [PIN,TAUN] = ANGDEPFUN_MIE(THETA,N) calculates the angle dependent
%   functions pi_n and tau_n (PIN, TAUN) for all angles in array THETA
%   and integer orders in array N.
%
%   Copyright 2012 Jan Schäfer, Institut für Lasertechnologien (ILM)
%   Author: Jan Schäfer (jan.schaefer@ilm.uni-ulm.de)
%   Organization: Institut für Lasertechnologien in der Medizin und
%       Meßtechnik an der Universität Ulm (http://www.ilm-ulm.de)

pin = zeros(numel(n), numel(theta));
taun = zeros(numel(n), numel(theta));
mu = cos(theta);
pin(1,:) = 1;
pin(2,:) = 3*mu;
taun(1,:) = mu;
taun(2,:) = 6*mu.^2 - 3;
for in=3:n(end)
    pin(in,:) = (2*in - 1)/(in - 1)*mu.*pin(in-1,:) ...
        - in/(in - 1).*pin(in-2,:);
    taun(in,:) = in*mu.*pin(in,:) - (in + 1)*pin(in-1,:);
end %for in=3:n(end)

end


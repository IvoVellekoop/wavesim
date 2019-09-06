function Ediff = f_mix(Ediff,Eprop,gamma)
% function applies medium potential gamma to propagated field Eprop and 
% previous field Ediff and combines field. 
Ediff = (1-gamma) .* Ediff - 1.0i * gamma.^2 .* Eprop;
end


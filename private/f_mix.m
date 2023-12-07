function Ediff = f_mix(Ediff,Eprop,gamma)
    % a mixing step:        E_diff  => (1-α γ) E_diff - i α γ² Eprop
    alpha = 0.7;
    Ediff = (1-alpha * gamma) .* Ediff - 1.0i * alpha * gamma.^2 .* Eprop;
end


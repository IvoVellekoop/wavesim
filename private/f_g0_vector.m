function [fEx, fEy, fEz] = f_g0_vector(fEx, fEy, fEz, px, py, pz, k0e)
% helper function for WaveSimVector. Calculates the dyadic Green's
% function for vector field propagation.
    div = (px .* fEx + py .* fEy + pz .* fEz) / k0e; %divergence term
    g0_k = 1 ./ (px.^2+py.^2+pz.^2-k0e);
    fEx = g0_k .* (fEx - px .* div);
    fEy = g0_k .* (fEy - py .* div);
    fEz = g0_k .* (fEz - pz .* div);
end

function E = f_iwiggle(E, g1, g2, g3)
% f_iwiggle Multiply E by the conjugate of the three gradient vectors g1, g2, g3. 
% method is used by both x-space and k-space wiggle algorithm
%todo: somehow this step is really slow even on the gpu!
    E = E .* conj(g1) .* conj(g2) .* conj(g3);
end


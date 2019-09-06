function E = f_wiggle(E, g1, g2, g3)
    % f_wiggle Multiply E by the three gradient vectors g1, g2, g3. method is
    % used by both x-space and k-space wiggle algorithm
    %todo: somehow this step is really slow even on the gpu!
    E = E .* g1 .* g2 .* g3;
end


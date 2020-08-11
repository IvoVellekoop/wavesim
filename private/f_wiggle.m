function E = f_wiggle(E, g1, g2, g3)
    % f_wiggle Multiply E by the three gradient vectors g1, g2, g3. method 
    % is used to perform the anti-cyclic convolution
    %todo: somehow this step is really slow even on the gpu!
    E = E .* g1 .* g2 .* g3;
end


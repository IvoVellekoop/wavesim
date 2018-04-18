function E = f_wiggle(E, gx, gy, gz)
%f_wiggle Multiply E by the three gradient vectors gx, gy, gz
    E = E .* gx .* gy .* gz;
end


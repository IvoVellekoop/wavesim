function E = f_wiggle(E, gx, gy, gz)
%f_wiggle Multiply E by the three gradient vectors gx, gy, gz
%todo: somehow this step is really slow even on the gpu!
%todo: ideally this could be implemented in the fft itself
    E = E .* gx .* gy .* gz;
end


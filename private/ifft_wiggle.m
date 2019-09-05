function E = ifft_wiggle(E, wig, gpu_enabled)
% Modified inverse Fourier transform: additionally applies wiggle
% phase ramps in real-space and k-space in wiggle descriptor wig.
% Todo: move to separate wiggle class?
if gpu_enabled
    E = arrayfun(@f_iwiggle, E, wig.gpx, wig.gpy, wig.gpz);
    E = ifftn(E);
    E = arrayfun(@f_iwiggle, E, wig.gx, wig.gy, wig.gz);
else
    E = f_iwiggle(E, wig.gpx, wig.gpy, wig.gpz);
    E = ifftn(E);
    E = f_iwiggle(E, wig.gx, wig.gy, wig.gz);
end
end
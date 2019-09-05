function E = fft_wiggle(E, wig, gpu_enabled)
% Modified Fourier transform: additionally applies wiggle phase
% ramps in real-space and k-space in wiggle descriptor wig.
% Todo: move to separate wiggle class?
if gpu_enabled
    E = arrayfun(@f_wiggle, E, wig.gx, wig.gy, wig.gz);
    E = fftn(E);
    %                 E = arrayfun(@f_wiggle, E, wig.gpx, wig.gpy, wig.gpz); % this line is not needed for the current implementation
else
    E = f_wiggle(E, wig.gx, wig.gy, wig.gz);
    E = fftn(E);
    %                 E = f_wiggle(E, wig.gpx, wig.gpy, wig.gpz); % this line is not needed for the current implementation
end
end
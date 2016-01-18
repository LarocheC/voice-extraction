%% Generate probe signal
probe = generate_probe(6554);
x = probe(1:32768).';

%% Setup filter banks
clear opts;
opts{1}.time.T = 256;
opts{1}.time.nFilters_per_octave = 8;
opts{1}.time.has_duals = true;
opts{2}.time = struct();
opts{2}.time.nFilters_per_octave = 1;
opts{2}.time.is_chunked = false;
archs = sc_setup(opts);

%% Compute scattering transform
[S, U, Y] = sc_propagate(x, archs);
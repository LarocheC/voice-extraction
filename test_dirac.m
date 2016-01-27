x = zeros(32768, 1);
x(16384) = 1;

%% Setup filter banks
clear opts;
phi_log2_oversampling = 2;
phi_oversampling = pow2(phi_log2_oversampling);
nmf_log2_oversampling = 2;

opts{1}.time.T = 1024;
opts{1}.time.nFilters_per_octave = 8;
opts{1}.time.has_duals = true;
opts{1}.time.is_chunked = false;
opts{1}.time.is_phi_gaussian = true;
opts{1}.time.size = 32768;
opts{2}.time.has_duals = true;
opts{2}.time.nFilters_per_octave = 2;
opts{2}.time.cutoff_in_dB = 2;
opts{2}.time.max_scale = length(x);
opts{2}.time.is_phi_gaussian = true;
opts{2}.time.S_log2_oversampling = phi_log2_oversampling;
archs = sc_setup(opts);

%% Compute scattering transform
[S, U, Y] = sc_propagate(x, archs);

x_test = mask_scattering(S{1+1}, S{1+2}, U, Y, archs);

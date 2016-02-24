%% Setup filter banks
addpath(genpath('..'));
clear opts;
phi_log2_oversampling = 6;
phi_oversampling = pow2(phi_log2_oversampling);
nmf_log2_oversampling = 3;

N = 32768;
T = N;
nFilters_per_octave = 8;

archs = make_archs(N, T, nFilters_per_octave, phi_log2_oversampling);
archs{1}.banks{1}.behavior.U.log2_oversampling = Inf;

%% Generate probe signal
[probe, noise, harmonic] = generate_probe(6554);
x = probe(1:32768).';
x = x - mean(x);

% Let's try a wavelet transform with full oversampling
arch = archs{1};
U0 = initialize_U(x, archs{1}.banks{1});

Y1 = U_to_Y(U0, archs{1}.banks);

%% Take the Fourier transform of each wavelet subband
Y1{end} = perform_ft(Y1{end}, archs{1}.banks{1}.behavior.key);

%% Initialize reconstruction to zero
rec_Y0 = Y1{1};
rec_Y0.data = zeros(N, 1);
rec_Y0.data_ft = zeros(N, 1);
%%

rec_Y0 = ...
    dual_scatter_dY(Y1{end}, arch.banks{1}, rec_Y0);

%%
rec_x = real(rec_Y0.data);
%% Setup filter banks
addpath(genpath('..'));
clear opts;
phi_log2_oversampling = 6;
phi_oversampling = pow2(phi_log2_oversampling);
nmf_log2_oversampling = 3;

N = 32768;
T = N;
nFilters_per_octave = 8;

opts{1}.time.T = T;
opts{1}.time.nFilters_per_octave = nFilters_per_octave;
opts{1}.time.has_duals = true;
opts{1}.time.is_chunked = false;
opts{1}.time.size = N;
opts{1}.time.S_log2_oversampling = phi_log2_oversampling;
opts{1}.time.U_log2_oversampling = Inf;
opts{1}.time.is_U_blurred = true;
opts{1}.time.trim_threshold = 0.0;

archs = sc_setup(opts);

%% Generate probe signal
% Real Signal
[noise, fs] = ...
    audioread('MusicDelta_FunkJazz_STEMS/MusicDelta_FunkJazz_STEM_01.wav');

[harmonic1, fs] = ...
    audioread('MusicDelta_FunkJazz_STEMS/MusicDelta_FunkJazz_STEM_02.wav');
[harmonic2, fs] = ...
    audioread('MusicDelta_FunkJazz_STEMS/MusicDelta_FunkJazz_STEM_03.wav');

noise = mean(noise,2);
noise = noise(1:N);

harmonic = mean(harmonic1+harmonic2,2);
harmonic = harmonic(1:N);

x = noise + harmonic;
x = x - mean(x);
%% Wavelet transform, forward and backward
arch = archs{1};
U0 = initialize_U(x, archs{1}.banks{1});

% Forward wavelet transform
Y1 = U_to_Y(U0, archs{1}.banks);

% Inverse wavelet transform
rec_Y0 = Y1{1}{1};
Y1{end} = perform_ft(Y1{end}, archs{1}.banks{1}.behavior.key);
rec_Y0_from_phi = dual_blur_dY(Y1{end}{2}, arch.banks{1});
rec_Y0.data_ft = 0.32 * rec_Y0_from_phi.data_ft;
rec_Y0 = ...
    dual_scatter_dY(Y1{end}{1}, arch.banks{1}, rec_Y0);

% Measure SNR
rec_x = real(rec_Y0.data);

plot(rec_x - x);
bss_eval_sources(rec_x',x')
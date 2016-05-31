%% Setup filter banks
clear opts;
phi_log2_oversampling = 9;
phi_oversampling = pow2(phi_log2_oversampling);
nmf_log2_oversampling = 2;

N = 2^17;
T = 2048;
nFilters_per_octave = 8;

archs = make_archs(N, T, nFilters_per_octave, phi_log2_oversampling);

%% Generate probe signal

% Real Signal
drum_path = 'MusicDelta_FunkJazz_STEMS/MusicDelta_FunkJazz_STEM_01.wav';
harmonic1_path = 'MusicDelta_FunkJazz_STEMS/MusicDelta_FunkJazz_STEM_02.wav';
harmonic2_path = 'MusicDelta_FunkJazz_STEMS/MusicDelta_FunkJazz_STEM_03.wav';

[drum_signal, sr] = audioread(drum_path);
harmonic1_signal = audioread(harmonic1_path);
harmonic2_signal = audioread(harmonic2_path);

drum_signal = mean(drum_signal, 2);
drum_signal = drum_signal(1:N);

harmonic_signal = mean(harmonic1_signal + harmonic2_signal,2);
harmonic_signal = harmonic_signal(1:N);

x = drum_signal + harmonic_signal;

[S, U, Y] = scattering(x, archs, phi_log2_oversampling, nmf_log2_oversampling);
Smat_percOracle = CreateDictionaryScattering( ...
    drum_signal, archs, phi_log2_oversampling, nmf_log2_oversampling,'all');
Smat_harmoOracle = CreateDictionaryScattering( ...
    harmonic_signal, archs, phi_log2_oversampling, nmf_log2_oversampling,'all');
load('xtrain');
xtrain = xtrain(1:N);
Smat_Train = CreateDictionaryScattering( ...
    xtrain', archs, phi_log2_oversampling, nmf_log2_oversampling,'all');


archs = downsample_archs(archs, nmf_log2_oversampling);
[Smat_1, Smat_2, refs_1, refs_2] = format_scattering(S);


subplot(121); imagesc((Smat_1));
subplot(122); imagesc((Smat_2));
%set(gcf, 'WindowStyle', 'docked');

%% % Launch SPNMF With a noise dictionary
% Create parameters
[F1, T] = size(Smat_1);
F2 = size(Smat_2, 1);
F = F1 + F2;

% W_perc = Smat_percOracle.^2;
W_perc = Smat_Train.^2;
H_perc = rand(size(W_perc,2), T);
max_iter = 200;
Smat = cat(1, Smat_1, Smat_2);
W_harmo = rand(F,50);

% SPNMF with no training permute the two layer in the synthetic case. 

[W_harmo, H_perc, e] = ...
    SPNMF_IS_W2TRAIN(Smat.^2, W_harmo, W_perc, H_perc, max_iter);

subplot(111);
plot(e);



%% Begin reconstruction 

% Smat_perc = Smat_percOracle; % case oracle
% Smat_harmo = Smat_harmoOracle; % case oracle

Smat_harmo = sqrt(W_harmo*W_harmo'*Smat);
Smat_perc = sqrt(W_perc * H_perc);
subplot(121); imagesc((Smat_harmo));
subplot(122); imagesc(((Smat_perc)));

Smat1_perc = Smat_perc(1:F1, :);
Smat2_perc = Smat_perc((F1+1):end, :);

Smat1_harmo = Smat_harmo(1:F1, :);
Smat2_harmo = Smat_harmo((F1+1):end, :);

% Fill in scattering network corresponding to harmonic part
S1_harmo = S{1+1};
for ref_1_index = 1:length(refs_1)
    ref_1 = refs_1(:, ref_1_index);
    row = Smat1_harmo(ref_1_index, :);
    S1_harmo.data = subsasgn(S1_harmo.data, ref_1, row);
end

S2_harmo = S{1+2};
for ref_2_index = 1:length(refs_2)
    ref_2 = refs_2(:, ref_2_index);
    row = Smat2_harmo(ref_2_index, :);
    S2_harmo.data = subsasgn(S2_harmo.data, ref_2, row);
end

% Fill in scattering network corresponding to percussive part
S1_perc = S{1+1};
for ref_1_index = 1:length(refs_1)
    ref_1 = refs_1(:, ref_1_index);
    row = Smat1_perc(ref_1_index, :);
    S1_perc.data = subsasgn(S1_perc.data, ref_1, row);
end

S2_perc = S{1+2};
for ref_2_index = 1:length(refs_2)
    ref_2 = refs_2(:, ref_2_index);
    row = Smat2_perc(ref_2_index, :);
    S2_perc.data = subsasgn(S2_perc.data, ref_2, row);
end

%%
x_perc = mask_scattering(S1_perc, S2_perc, U, Y, archs);
x_harmo = mask_scattering(S1_harmo, S2_harmo, U, Y, archs);


%% Compute BSS score

se = [x_perc' ; x_harmo'];
s = [drum_signal' ; x'-drum_signal'];

[SDR,SIR,SAR,perm]=bss_eval_sources(se,s);


%%
% subplot(131); imagesc((W_perc*H_perc));
% subplot(132); imagesc(Smat);


%% STFT comparison 

% parameters
M = 1024;
a = M/2;
g = gabwin({'tight', 'hann'}, a, M);
Ls = length(x);
X = dgtreal(x,g,a,M);
Xnoise = dgtreal(drum_signal,g,a,M);
Xharmo = dgtreal(harmonic_signal,g,a,M);


[F, T] = size(X);
Wini = rand(F,50);

Wtrain = abs(dgtreal(xtrain,g,a,M));
scale = sqrt(sum(Wtrain.^2,1)); 
Wtrain = Wtrain .* repmat((scale+eps).^-1,F,1);
WP =  Wtrain.^2;
HPini = rand(size(WP,2),T);


max_iter = 200;

% SPNMF with no training permute the two layer in the synthetic case. 

[W, HP, estft] = ...
    SPNMF_IS_W2TRAIN(abs(X).^2, Wini, WP, HPini, max_iter);


subplot(121); imagesc(10*log10(W*W'*abs(X)),[-50 10]);
subplot(122); imagesc(10*log10(WP * HP),[-50 10])


% Wiener Filtering and reconstruction of the signal
V_harmo = (W*W'*abs(X)).^2;
V_perc = (WP*HP).^2;

% V_harmo = abs(Xharmo).^2; %oracle case
% V_perc = abs(Xnoise).^2; %oracle case

V = V_harmo + V_perc;

DrumEst = idgtreal((V_perc./(V+eps)).*(X),g,a,M,Ls);

HarmoEst = idgtreal((V_harmo./(V+eps)).*(X),g,a,M,Ls);

[SDRstft,SIRstft,SARstft,perm]=bss_eval_sources([DrumEst' ; HarmoEst'], s);

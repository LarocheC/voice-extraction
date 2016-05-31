% %% Setup filter banks
% clear opts;
phi_log2_oversampling = 9;
phi_oversampling = pow2(phi_log2_oversampling);
nmf_log2_oversampling = 2;

load('hits')

%% Real Signal
drum_path = 'MusicDelta_FunkJazz_STEMS/MusicDelta_FunkJazz_STEM_01.wav';
harmonic1_path = 'MusicDelta_FunkJazz_STEMS/MusicDelta_FunkJazz_STEM_03.wav';
harmonic2_path = 'MusicDelta_FunkJazz_STEMS/MusicDelta_FunkJazz_STEM_04.wav';

[drum_signal, sr] = audioread(drum_path);
harmonic1_signal = audioread(harmonic1_path);
harmonic2_signal = audioread(harmonic2_path);

drum_signal = mean(drum_signal, 2);
%drum_signal = drum_signal(1:N);

harmonic_signal = mean(harmonic1_signal + harmonic2_signal,2);
%harmonic_signal = harmonic_signal(1:N);

x = drum_signal + harmonic_signal;

N = 2^20;
T = 2^17;
nFilters_per_octave = 8;

archs = make_archs(N, T, nFilters_per_octave, phi_log2_oversampling);


[S, U, Y] = scattering(x(1:N), archs, phi_log2_oversampling, nmf_log2_oversampling);
[Smat_1x, Smat_2x, refs_1, refs_2] = format_scattering(S);

drummer_index = 2;
hit_index = 1;
hit_file_index = [1 4 5 6 7 8];

archs_Dictionary = make_archs(T, T, nFilters_per_octave, phi_log2_oversampling);

for index=1:6
    [S, U, Y] = scattering(hits{drummer_index}{hit_file_index(index)}{hit_index},...
        archs_Dictionary, phi_log2_oversampling, nmf_log2_oversampling);
    [Smat_1, Smat_2, refs_1, refs_2] = format_scattering(S);
    dictionary(index).Smat = cat(1, Smat_1, Smat_2);
    dictionary(index).U = U;
    dictionary(index).Y = Y;
    dictionary(index).S = S;

    if index==1
        Smat_Dict = dictionary(1).Smat;
    else
        Smat_Dict = [Smat_Dict dictionary(index).Smat];
    end

end


%% % Launch SPNMF With a noise dictionary
% Create parameters
[F1, T] = size(Smat_1x);
F2 = size(Smat_2x, 1);
F = F1 + F2;

% W_perc = Smat_percOracle.^2;
W_perc = Smat_Dict.^2;
H_perc = rand(size(W_perc,2), T);
max_iter = 200;
Smat = cat(1, Smat_1x, Smat_2x);
W_harmo = rand(F,50);

% SPNMF with no training permute the two layer in the synthetic case. 

[W_harmo, H_perc, e] = ...
    SPNMF_IS_W2TRAIN(Smat.^2, W_harmo, W_perc, H_perc, max_iter);

subplot(111);
plot(e);

Smat_harmo = W_harmo*W_harmo'*Smat;
Smat_perc = W_perc * H_perc;
subplot(121); imagesc(db(Smat_harmo),[-120 60]);
subplot(122); imagesc(db(Smat_perc));






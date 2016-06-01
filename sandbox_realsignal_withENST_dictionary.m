% %% Setup filter banks
% clear opts;
phi_log2_oversampling = 9;
phi_oversampling = pow2(phi_log2_oversampling);
nmf_log2_oversampling = 2;

load('hits')

N = 2^20;

%% Real Signal
drum_path = 'MusicDelta_FunkJazz_STEMS/MusicDelta_FunkJazz_STEM_01.wav';
harmonic1_path = 'MusicDelta_FunkJazz_STEMS/MusicDelta_FunkJazz_STEM_03.wav';
harmonic2_path = 'MusicDelta_FunkJazz_STEMS/MusicDelta_FunkJazz_STEM_04.wav';

[drum_signal, sr] = audioread(drum_path);
harmonic1_signal = audioread(harmonic1_path);
harmonic2_signal = audioread(harmonic2_path);

drum_signal = mean(drum_signal, 2);
drum_signal = drum_signal(1:N);

harmonic_signal = mean(harmonic1_signal + harmonic2_signal,2);
harmonic_signal = harmonic_signal(1:N);

x = drum_signal + harmonic_signal;

T = 2^17;
nFilters_per_octave = 8;

archs = make_archs(N, T, nFilters_per_octave, phi_log2_oversampling);

[S, U, Y] = scattering(x, archs, phi_log2_oversampling, nmf_log2_oversampling);
[Smat_1x, Smat_2x, refs_1, refs_2] = format_scattering(S);

% [S, U, Y] = scattering(harmonic_signal(1:N), archs, phi_log2_oversampling, nmf_log2_oversampling);
% [Smat_1Harmo, Smat_2Harmo, refs_1, refs_2] = format_scattering(S);
% 
% [S, U, Y] = scattering(drum_signal(1:N), archs, phi_log2_oversampling, nmf_log2_oversampling);
% [Smat_1Percu, Smat_2Percu, refs_1, refs_2] = format_scattering(S);



drummer_index = 2;
hit_index = 1;
hit_file_index = [1 4 5 6 7 8];

archs_Dictionary = make_archs(T, T, nFilters_per_octave, phi_log2_oversampling);

for index=1:6
    [S1, U1, Y1] = scattering(hits{drummer_index}{hit_file_index(index)}{hit_index},...
        archs_Dictionary, phi_log2_oversampling, nmf_log2_oversampling);
    [Smat_1, Smat_2, refs_1, refs_2] = format_scattering(S1);
    dictionary(index).Smat = cat(1, Smat_1, Smat_2);
    dictionary(index).U = U1;
    dictionary(index).Y = Y1;
    dictionary(index).S = S1;

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

W_perc = Smat_Dict;
H_perc = rand(size(W_perc,2), T);
max_iter = 200;
Smat = cat(1, Smat_1x, Smat_2x);
% SmatHarmo = cat(1, Smat_1Harmo, Smat_2Harmo);
% SmatPercu = cat(1, Smat_1Percu, Smat_2Percu);
W_harmo = rand(F,50);

% SPNMF 

[W_harmo, H_perc, e] = ...
    SPNMF_KL_W2TRAIN(Smat, W_harmo, W_perc, H_perc, max_iter);

plot(e);

Smat_harmo = W_harmo*W_harmo'*Smat;
Smat_perc = W_perc * H_perc;

% subplot(221); imagesc(db(SmatHarmo),[-120 60]);
% subplot(222); imagesc(db(SmatPercu),[-120 60]);
% subplot(223); imagesc(db(Smat_harmo),[-120 60]);
% subplot(224); imagesc(db(Smat_perc),[-120 60]);


%% Begin reconstruction 


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



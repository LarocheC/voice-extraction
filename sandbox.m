%% Generate probe signal
[probe, noise] = generate_probe(6554);
x = probe(1:32768).';
noise = noise(1:32768).';

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


% Generate references
spatial_subscripts = 1;
refs_1 = generate_refs(S{1+1}.data, spatial_subscripts, S{1+1}.ranges{1+0});
refs_2 = generate_refs(S{1+2}.data, spatial_subscripts, S{1+2}.ranges{1+0});

% Convert scattering nodes to matrices
nFrames = size(S{1+1}.data, 1);
nRefs_1 = size(refs_1, 2);
Smat_1 = zeros(nRefs_1, nFrames);
for ref1 = 1:nRefs_1
    Smat_1(ref1, :) = subsref(S{1+1}.data, refs_1(ref1));
end
nRefs_2 = size(refs_2, 2);
Smat_2 = zeros(nRefs_2, nFrames);
for ref2 = 1:nRefs_2
    Smat_2(ref2, :) = subsref(S{1+2}.data, refs_2(:, ref2));
end

% Downsample manually
nmf_oversampling = pow2(nmf_log2_oversampling);
downsampling = phi_oversampling / nmf_oversampling;
Smat_1 = Smat_1(:, 1:downsampling:end);
Smat_2 = Smat_2(:, 1:downsampling:end);

% Update the toolbox behavior so that it adapts to resampled inputs
archs{2}.banks{1}.behavior.S.log2_oversampling = nmf_log2_oversampling;
archs{2}.invariants{1}.behavior.S.log2_oversampling = nmf_log2_oversampling;
S{1+1}.ranges{1+0}(2,1) = S{1+1}.ranges{1+0}(2,1) * downsampling;

% Ensure non-negativity
Smat_1(Smat_1<0)=0;
Smat_2(Smat_2<0)=0;

subplot(121); imagesc((Smat_1));
subplot(122); imagesc((Smat_2));
%set(gcf, 'WindowStyle', 'docked');

%% % Launch SPNMF With a noise dictionary
% Create parameters
[F1, T] = size(Smat_1);
F2 = size(Smat_2, 1);
F = F1 + F2;

W_harmo = rand(F, 4);
W_perc = rand(F, 2);
H_perc = rand(2, T);
max_iter = 200;

Smat = cat(1, Smat_1, Smat_2);

% SPNMF with no training permute the two layer in the synthetic case. 

[W_harmo, W_perc, H_perc, e] = ...
    SPNMF_KL(Smat, W_harmo, W_perc, H_perc, max_iter);

subplot(111);
plot(e);

subplot(121); imagesc(W_harmo*W_harmo'*Smat,[0 35]);
subplot(122); imagesc(W_perc * H_perc)

Smat_perc = W_perc * H_perc;
Smat1_perc = Smat_perc(1:F1, :);
Smat2_perc = Smat_perc((F1+1):end, :);

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
Y2_perc = dS_backto_dY(S1_perc, archs{2});
Y2_perc = copy_metadata(Y{1+1}{1}, Y2_perc);
U1fromS1_perc = dY_backto_dU(Y2_perc);
U1fromS1_perc = perform_ft(U1fromS1_perc, archs{1}.banks{1}.behavior.key);

Y3_perc = dS_backto_dY(S2_perc, archs{3});
Y3_perc = copy_metadata(Y{1+2}{1}, Y3_perc);
U2_perc = dY_backto_dU(Y3_perc);
Y2_perc = Y{2}{end};
for lambda2 = 1:length(U2_perc.data)
    Y2_perc.data{lambda2} = U2_perc.data{lambda2} .* ...
        Y{2}{end}.data{lambda2} ./ abs(Y{2}{end}.data{lambda2});
end

Y1_perc = dual_scatter_dY(Y2_perc, archs{2}.banks{1}, U1fromS1_perc);
U1_perc = dY_backto_dU(Y1_perc);
for lambda1 = 1:length(U1fromS1_perc.data)
    Y1_perc.data{lambda1} = U1_perc.data{lambda1} .* ...
        Y{1}{end}.data{lambda1} ./ abs(Y{1}{end}.data{lambda1});
end
%
Y0_perc = Y{1+0}{1};
Y0_perc.data = complex(zeros(size(Y0_perc.data)));
Y0_perc.data_ft = complex(zeros(size(Y0_perc.data_ft)));
Y0_perc = dual_scatter_dY(Y1_perc, archs{1}.banks{1}, Y0_perc);
x_perc = real(Y0_perc.data);


%% Compute BSS score

se = [x_perc' ; x'-x_perc'];
s = [noise' ; x'-noise'];

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


[F, T] = size(X);
Wini = rand(F, 4);
WPini = rand(F, 2);
HPini = rand(2, T);
max_iter = 200;

% SPNMF with no training permute the two layer in the synthetic case. 

[W, WP, HP, e] = ...
    SPNMF_KL(abs(X).^2, Wini, WPini, HPini, max_iter);


subplot(121); imagesc(10*log10(W*W'*abs(X)),[-50 10]);
subplot(122); imagesc(10*log10(WP * HP),[-50 10])


% Wiener Filtering and reconstruction of the signal

V = (W*W'*abs(X)).^2 + (WP*HP).^2 ;

DrumEst = idgtreal(((WP*HP).^2./(V+eps)).*(X),g,a,M,Ls);

HarmoEst = idgtreal((((W*W'*abs(X)).^2)./(V+eps)).*(X),g,a,M,Ls);

[SDRstft,SIRstft,SARstft,perm]=bss_eval_sources([DrumEst' ; HarmoEst'], s);



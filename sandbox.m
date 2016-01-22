close all


%% Generate probe signal
probe = generate_probe(6554);
x = probe(1:32768).';


%% Setup filter banks
clear opts;
opts{1}.time.T = 1024;
opts{1}.time.nFilters_per_octave = 8;
opts{1}.time.has_duals = true;
opts{1}.time.is_chunked = false;
opts{1}.time.size = 32768;
opts{1}.is_phi_gaussian = true;
opts{2}.is_phi_gaussian = true;
opts{2}.time = struct();
opts{2}.time.nFilters_per_octave = 1;
archs = sc_setup(opts);

%% Compute scattering transform
[S, U, Y] = sc_propagate(x, archs);

%% Generate references
spatial_subscripts = 1;
refs_1 = generate_refs(S{1+1}.data, spatial_subscripts, S{1+1}.ranges{1+0});
refs_2 = generate_refs(S{1+2}.data, spatial_subscripts, S{1+2}.ranges{1+0});

%% Convert scattering nodes to matrices
nFrames = size(S{1+1}.data, 1);
nRefs_1 = size(refs_1, 2);
Smat_1 = zeros(nFrames, nRefs_1);
for ref1 = 1:nRefs_1
    Smat_1(:, ref1) = subsref(S{1+1}.data, refs_1(:, ref1));
end
nRefs_2 = size(refs_2, 2);
Smat_2 = zeros(nFrames, nRefs_2);
for ref2 = 1:nRefs_2
    Smat_2(:, ref2) = subsref(S{1+2}.data, refs_2(:, ref2));
end

Smat_1(Smat_1<0)=0;
Smat_2(Smat_2<0)=0;

Smat_1 = Smat_1./(norm(Smat_1,2)+eps);
Smat_2(Smat_2<0)=0;

% subplot(121); imagesc((Smat_1'));
% subplot(122); imagesc((Smat_2'));


%% Create parameters
[F1,T1]=size(Smat_1');
[F2,T2]=size(Smat_2');


W = rand(F1,4);
W2 = rand(F1,2);
H2 = rand(2,T1);
max_iter = 1000;

X = Smat_1';

%% Lunch SPNMF With a noise dictionary

[ W, H2, e] = SPNMF_KL_W2TRAIN( X , W , W2, H2 , max_iter);

plot(e);
 
figure
subplot(131); imagesc((W*W'*Smat_1'));
subplot(132); imagesc((W2*H2));
subplot(133); imagesc((W*W'*Smat_1'+W2*H2));



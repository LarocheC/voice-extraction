%% Generate probe signal
probe = generate_probe(6554);
x = probe(1:32768).';

%% Setup filter banks
clear opts;
opts{1}.time.T = 8192;
opts{1}.time.nFilters_per_octave = 8;
opts{1}.time.has_duals = true;
opts{1}.time.is_chunked = false;
opts{1}.time.is_phi_gaussian = true;
opts{1}.time.size = 32768;
opts{2}.time.nFilters_per_octave = 1;
opts{2}.time.is_phi_gaussian = true;
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
% For Vincent
% enst_drums_path = fullfile('~', 'datasets', 'ENST-drums-public');

% For Clement
enst_drums_path = 'C:\Users\laroche\Desktop\ENST-drums-public';

%% Segment hits
hits = segment_hits(enst_drums_path);

%%
nmf_log2_oversampling = 3;
phi_log2_oversampling = 6;

N = 131072;
T = 1024;
nFilters_per_octave = 8;

archs = make_archs(N, T, nFilters_per_octave, phi_log2_oversampling);

%%
nHits = length(hits);
atoms = cell(1, nHits);
planck_window = planck_taper(N);

for hit_index = 1:nHits
    disp(hit_index);
    hit = hits{hit_index};
    windowed_hit = hit .* planck_window;
    S = sc_propagate(windowed_hit, archs);
    S1 = S{1+1}.data.';
    S2 = [S{1+2}.data{:}].';
    spectral_flux = sum(S1);
    modulation_flux = sum(S2);
    [~, onset] = max(spectral_flux + modulation_flux);
    atoms{hit_index} = cat(1, S1(:, onset), S2(:, onset));
end
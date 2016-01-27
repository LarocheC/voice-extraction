function [archs, opts] = make_archs(N, T, nFilters_per_octave, phi_log2_oversampling)
    %phi_log2_oversampling = 6;%log2(N/T); % check me
    opts{1}.time.T = T;
    opts{1}.time.nFilters_per_octave = nFilters_per_octave;
    opts{1}.time.has_duals = true;
    opts{1}.time.is_chunked = false;
    opts{1}.time.is_phi_gaussian = true;
    opts{1}.time.size = N;
    opts{1}.time.S_log2_oversampling = phi_log2_oversampling;
    opts{2}.time.has_duals = true;
    opts{2}.time.nFilters_per_octave = 2;
    opts{2}.time.cutoff_in_dB = 2;
    opts{2}.time.max_scale = N;
    opts{2}.time.is_phi_gaussian = true;
    opts{2}.time.S_log2_oversampling = phi_log2_oversampling;
    archs = sc_setup(opts);
end
% For Vincent
enst_drums_path = fullfile('~', 'datasets', 'ENST-drums-public');

%%
drummer_str = 'drummer_1';
mix_str = 'wet_mix';
mix_path = fullfile(enst_drums_path, drummer_str, mix_str);

mix_path = fullfile(enst_drums_path, drummer_str, 'audio', mix_str);
dir_content = dir(mix_path);
all_file_names = {dir_content(~[dir_content.isdir]).name};
hit_file_names = ...
    all_file_names(cellfun( @(x) strcmp(x(5:8), 'hits'), all_file_names));
nHit_files = length(hit_file_names);

%for hit_file_index = 1:nHit_files
hit_file_index = 1;
    hit_file_name = hit_file_names{hit_file_index};
    hit_file_path = fullfile(mix_path, hit_file_name);
    x_indices = findstr(hit_file_path, 'x');
    last_x_index = x_indices(end);
    nHits = str2num(hit_file_path((last_x_index+1):(end-4)))
    stereo_waveform = audioread(hit_file_path);
    mono_waveform = mean(stereo_waveform, 2);
    [peak_values, peak_locations] = findpeaks(abs(mono_waveform), ...
        'NPeaks', nHits, 'sort', 'descend', 'MinPeakDistance', 10000)
%end
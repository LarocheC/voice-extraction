% For Vincent
enst_drums_path = fullfile('~', 'datasets', 'ENST-drums-public');

%%
nDrummers = 3;
hits = cell(1, nDrummers);

for drummer_index = 1:nDrummers
    drummer_str = ['drummer_', num2str(drummer_index)];
    mix_str = 'wet_mix';
    mix_path = fullfile(enst_drums_path, drummer_str, 'audio', mix_str);
    
    dir_content = dir(mix_path);
    all_file_names = {dir_content(~[dir_content.isdir]).name};
    hit_file_names = ...
        all_file_names(cellfun( @(x) strcmp(x(5:8), 'hits'), all_file_names));
    nHit_files = length(hit_file_names);
    
    nfft = 1024;
    hit_length = 131072;
    hits{drummer_index} = cell(1, nHit_files);
    
    for hit_file_index = 1:nHit_files
        hit_file_name = hit_file_names{hit_file_index};
        hit_file_path = fullfile(mix_path, hit_file_name);
        x_indices = strfind(hit_file_path, 'x');
        last_x_index = x_indices(end);
        nHits = str2double(hit_file_path((last_x_index+1):(end-4)));
        [stereo_waveform, sample_rate] = audioread(hit_file_path);
        mono_waveform = mean(stereo_waveform, 2);
        odf = onset_detection_function(mono_waveform, sample_rate, nfft);
        [~, hit_locations] = ...
            findpeaks(odf, ...
            'NPeaks', nHits, ...
            'sort', 'descend', ...
            'MinPeakDistance', 44100 / nfft);
        nHits = length(hit_locations);
        hit_locations = hit_locations * nfft;
        hit_locations = sort(hit_locations, 'ascend');
        hits{drummer_index}{hit_file_index} = cell(1, nHits);
        for hit_index = 1:nHits
            hit_location = hit_locations(hit_index);
            hit_start = hit_location -  hit_length / 4;
            if hit_index > 1
                hit_start = max(hit_start, hit_locations(hit_index-1));
            else
                hit_start = max(hit_start, 1);
            end
            hit_stop = hit_start + hit_length - 1;
            if hit_index < nHits
                hit_stop = min(hit_stop, ...
                    hit_locations(hit_index+1) - hit_length/4);
            else
                hit_stop = min(hit_stop, length(mono_waveform));
            end
            hit_waveform = mono_waveform(hit_start:hit_stop);
            if hit_start == 1
                hit_waveform = cat(1, zeros(hit_length/4 - 1, 1), ...
                    hit_waveform);
            end
            if length(hit_waveform) < hit_length
                hit_waveform = cat(1, hit_waveform, ...
                    zeros(hit_length - length(hit_waveform), 1));
            end
            hits{drummer_index}{hit_file_index}{hit_index} = hit_waveform;
        end
        hits{drummer_index} = [hits{drummer_index}{:}];
    end
    hits = [hits{:}];
end
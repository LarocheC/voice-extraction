function detection_function = onset_detection(waveform, sample_rate)
%% Compute spectrogram
nfft = 1024;
frame_duration = 10; % in ms
frame_length = floor(sample_rate * frame_duration/1000);
window = hamming(2*frame_length-1);
spectrum = spectrogram(waveform, window, frame_length, nfft, sample_rate);

%% Compute spectral flux
denominator = 10.0;
log_spectrum = log1p(abs(spectrum) / denominator);
spectral_flux = sum(abs(diff(log_spectrum, 2)), 1);

%% Low-pass filtering
lowpass_cutoff_length = 15; % in number of windows
hann_window = hanning(2*lowpass_cutoff_length);
lowpass_filter = hann_window(lowpass_cutoff_length+1:end);
detection_function = filter(lowpass_filter, 1, spectral_flux);
end


function [probe, noise] = generate_probe(fs)
if nargin<1
    fs = 4000;
end
%%
t = 0:1/fs:5; % 5 secs @ 4000Hz sample rate

T = 2.5;
T1 = 4;
f1 = 130.813; %C3
f3 = 220.000; %A3

amplitude_1 = max(t.^1.5 .* exp(-2 * t), 0);
amplitude_1 = amplitude_1 / max(amplitude_1);

amplitude_2 = max((t-2.2).^3 .* exp(-4.5 * (t-2.2)), 0);
amplitude_2 = amplitude_2 / max(amplitude_2);

harmonic_1 = 0.5 * amplitude_1 .* generate_harmonics(f3, t, 3); 
harmonic_2 = 0.5 * amplitude_2 .* generate_harmonics(f1, t, 3);

noise_1 = 2.0 * (heaviside(t-1) - heaviside(t-1.05)) .* (rand(1, length(t))-0.5);
noise_2 = 1.0 * (heaviside(t-2) - heaviside(t-2.1)) .* (rand(1, length(t))-0.5);
noise_3 = 0.5 * (heaviside(t-3) - heaviside(t-3.2)) .* (rand(1, length(t))-0.5);
noise_4 = 0.5 * (heaviside(t-4) - heaviside(t-4.4)) .* (rand(1, length(t))-0.5);

noise = noise_1 + noise_2 + noise_3 + noise_4;

harmochirp =       sin(2*pi*1*f1*(1+0.1*sin(2*t)) .* t) + ...
             0.3 * sin(2*pi*2*f1*(1+0.1*sin(2*t)) .* t) + ...
             0.1 * sin(2*pi*3*f1*(1+0.1*sin(2*t)) .* t);


chirp_1 = 0 * (heaviside(t-0.5) - heaviside(t-1.5)) .* harmochirp;
chirp_2 = 0 * (heaviside(t-2.5) - heaviside(t-3.5)) .* harmochirp;

probe = harmonic_1 + harmonic_2 + ...
        noise_1 + noise_2 + noise_3 + noise_4 + ...
        chirp_1 + chirp_2;
    plot(t,probe)
    soundsc(probe, fs)
end


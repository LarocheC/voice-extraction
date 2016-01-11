clear all 
close all

%% Signal de test multi couches

fs=4000;
t = 0:1/fs:5; % 5 secs @ 4000Hz sample rate

T = 2.5; T1 =4;
f1 = 131 ; f3 = 262 ;

x1 = 0.5*synth(f1,t,3).*(1-heaviside(t-3));
x2 = 0.5*(heaviside(t-2)).* synth(f3,t,3); 
x3 = 2*(heaviside(t-1)-heaviside(t-1.1)+heaviside(t-2)-heaviside(t-2.1))...
    .* rand(1,length(t)); %synth(f2,t,5);
x4 = 0.5*(heaviside(t-3)-heaviside(t-3.3)+heaviside(t-4)-heaviside(t-4.3)) .* rand(1,length(t));
x4 = (heaviside(t-0.5)-heaviside(t-1.5)+ heaviside(t-2.5)-heaviside(t-3.5)).*(sin(2* pi* f1* (1+0.1*sin(2*t)).* t) + 0.3*sin(2* pi* 2*f1* (1+0.1*sin(2*t)).* t)+ 0.1*sin(2* pi* 3*f1* (1+0.1*sin(2*t)).* t));


y = x4;               
b = fir1(34,0.48,'low',chebwin(35,30));    % FIR filter design
%freqz(b,1,512);                 % Frequency response of filter
x4 = filtfilt(b,1,y);       % Zero-phase digital filtering

y = x1;               
b = fir1(34,0.48,'low',chebwin(35,30));    % FIR filter design
%freqz(b,1,512);                 % Frequency response of filter
x1 = filtfilt(b,1,y);       % Zero-phase digital filtering

y = x2;               
b = fir1(34,0.48,'low',chebwin(35,30));    % FIR filter design
%freqz(b,1,512);                 % Frequency response of filter
x2 = filtfilt(b,1,y);       % Zero-phase digital filtering



x =x1+x2+x3+x4;
M = 256;
a = (M)/2;
g = gabwin({'tight', 'hann'},a,M);
X = dgtreal(x,g,a,M);

figure;
imagesc(10*log10(abs(X)),[-50 10]);


x = x.';
y = x(1000:7000);
figure
plot(y);
%%

% Setup
clear opts;
opts{1}.time.T = 256;
opts{1}.time.nFilters_per_octave = 8;
opts{1}.time.has_duals = true;
opts{2}.time = struct();
opts{2}.time.nFilters_per_octave = 1;
archs = sc_setup(opts);

% Propagate
[S, U, Y] = sc_propagate(x, archs);

% Plot S2 along with S1
subplot(211);
imagesc(S{1+1}.data.');

subplot(212);
imagesc([S{1+2}.data{:}].');

%% First-order display. Bank 1, U1, and S1

subplot(311);
display_bank(archs{1}.banks{1}, 'psis');

subplot(312);
U_unchunked = sc_unchunk(U);
display_scalogram(U_unchunked{1+1});

subplot(313);
imagesc(S{1+1}.data.');

%% Unchunk and format
S_unchunked = sc_unchunk(S);
S_matrix = sc_format(S_unchunked);


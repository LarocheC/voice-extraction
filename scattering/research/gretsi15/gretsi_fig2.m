% Reproduit la figure 2a de la soumission au GRETSI 2015
% Vincent Lostanlen, Stephane Mallat.
% "Transformee de scattering en spirale temps-chroma-octave"

%% Construction of scattering architectures, i.e. filterbanks and nonlinearities
% Number of samples
N = 16384;

% Options for scalogram
opts{1}.time.size = N;
opts{1}.time.nFilters_per_octave = 16;

% Options for scattering along time
opts{2}.time.handle = @gammatone_1d;
opts{2}.time.max_scale = Inf;
opts{2}.time.max_Q = 2;
opts{2}.time.U_log2_oversampling = 2;
opts{2}.time.gamma_bounds = [8 Inf];

% Options for scattering along chromas
opts{2}.gamma.invariance = 'bypassed';
opts{2}.gamma.U_log2_oversampling = Inf;

% Options for scattering along octaves
opts{2}.j.handle = @RLC_1d;
opts{2}.j.invariance = 'bypassed';
opts{2}.j.mother_xi = 0.5;

% Build scattering "architectures", i.e. filter banks and nonlinearities
archs = sc_setup(opts);

%% Computation of spiral scattering
% We start by computing an empty scalogram
signal = zeros(N,1);
U{1+0} = initialize_U(signal,archs{1}.banks{1});
Y{1} = U_to_Y(U{1+0},archs{1});
U{1+1} = Y_to_U(Y{1}{end},archs{1});

% We put a Dirac peak
U{1+1}.data{round(end/3)}(end/2) = 1;

% We compute the spiral scattering transform without modulus nonlinearity
Y{2} = U_to_Y(U{1+1},archs{2});

%% Visualization
% Settings for alpha, beta and gamma (see paper) are here
time_scale = 10;
chroma_scale = 1;
octave_scale = 2;
chroma_sign = 2;
octave_sign = 1;

% Coefficient extraction in a 2d
spiral_scattergram = ...
    Y{2}{1+3}{1,1,1}.data{time_scale}{chroma_scale,octave_scale};
sizes = size(spiral_scattergram);
scattergram_sizes = [sizes(1),sizes(2)*sizes(3),sizes(4),sizes(5)];
spiral_scattergram = reshape(spiral_scattergram,scattergram_sizes);
gamma_range = 1:(min(Inf,size(spiral_scattergram,2)));
spiral_scattergram = real(spiral_scattergram(:,gamma_range,:,:));
spiral_scattergram(abs(spiral_scattergram)<1e-4) = 0;
normalizer = max(abs(spiral_scattergram(:)));
spiral_scattergram = 32 + 32 * spiral_scattergram/normalizer;

% Display
colormap rev_hot;
image(spiral_scattergram(:,1:end/2,chroma_sign,octave_sign)');
axis off;

%% Export
%export_fig gretsi_fig2a.png -transparent

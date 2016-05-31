i = i+1;
plot(hits{i}); sound(hits{i}, sample_rate);


for i=1:4
    
    sel = [1 10 50 200];
    % parameters
    M = 1024;
    a = M/2;
    g = gabwin({'tight', 'hann'}, a, M);
    Ls = length(hits{sel(i)});
    X = dgtreal(hits{sel(i)},g,a,M);

    subplot(2,2,i)
    imagesc([0 Ls/sample_rate],[0 sample_rate/2],10*log10(abs(X)),[-50 10]); axis xy

end

%% The following indices should be discarded: 140, 142, 144, 420
discarded_indices = [140, 142, 144, 420];
hits(discarded_indices) = [];
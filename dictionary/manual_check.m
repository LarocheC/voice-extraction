i = i+1;
plot(hits{i}); sound(hits{i}, sample_rate);

%% The following indices should be discarded: 140, 142, 144, 420
discarded_indices = [140, 142, 144, 420];
hits(discarded_indices) = [];
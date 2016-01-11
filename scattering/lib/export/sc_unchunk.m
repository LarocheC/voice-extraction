function unchunked_S = sc_unchunk(S)
%% Initialization
nLayers = length(S);
unchunked_S = cell(1,nLayers);

%% Layer-wise unchunking
for layer_index = 0:nLayers-1
    unchunked_S{1+layer_index} = unchunk_layer(S{1+layer_index});
end
end
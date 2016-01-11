function layer_out = copy_metadata(layer_in,layer_out)
%% Deep map across cells
if iscell(layer_in)
    for cell_index = 1:numel(layer_in)
        if ~isempty(layer_out{cell_index})
            layer_out{cell_index} = ...
                copy_metadata(layer_in{cell_index},layer_out{cell_index});
        end
    end
else
    layer_out.keys = layer_in.keys;
    layer_out.ranges = layer_in.ranges;
    layer_out.variable_tree = layer_in.variable_tree;
end
end
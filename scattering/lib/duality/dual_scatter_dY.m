function previous_sub_dY = dual_scatter_dY(sub_dY,bank,previous_sub_dY)
%% Cell-wise map
if iscell(sub_dY)
    dual_scatter_handle = @(x,y) dual_scatter_dY(x,bank,y);
    previous_sub_dY = ...
        map_unary_inplace(dual_scatter_handle,sub_dY,previous_sub_dY);
    return
end

%% Variable loading
keys = sub_dY.keys;
ranges = sub_dY.ranges;
variable_tree = sub_dY.variable_tree;
[sibling,~,gamma_variable] = ...
    dual_get_relatives(bank.behavior.key,variable_tree);
variable = get_leaf(variable_tree,bank.behavior.key);

% Subscripts and colons are updated according to the network structure
bank.behavior.gamma_subscript = gamma_variable.subscripts;
bank.behavior.gamma_level = gamma_variable.level;
bank.behavior.subscripts = variable.subscripts;
bank.behavior.colons.subs = replicate_colon(length(keys{1+0}));

%% Dual scattering
if isempty(sibling)
    previous_sub_dY.data = ...
        dual_firstborn_scatter(sub_dY.data,bank,ranges, ...
        previous_sub_dY.data_ft,previous_sub_dY.ranges);
elseif sibling.nSiblings==0
    previous_sub_dY.data = ...
        dual_secondborn_scatter(sub_dY.data,bank,ranges,sibling, ...
        previous_sub_dY.data_ft,previous_sub_dY.ranges);
else
    previous_sub_dY.data = ...
        dual_sibling_scatter(sub_dY.data,bank,ranges,sibling, ...
        previous_sub_dY.data_ft,previous_sub_dY.ranges);
end
end
function nonlinearity = fill_nonlinearity(opt)
%%
found = false;
if isfield(opt,'nonlinearity')
    if isfield(opt.nonlinearity, 'denominator')
        opt.nonlinearity.name = ...
            default(opt.nonlinearity, 'name', 'uniform_log');
    end
    nonlinearity = opt.nonlinearity;
    if isfield(opt.nonlinearity,'name')
        switch opt.nonlinearity.name
            case 'modulus'
                nonlinearity.is_modulus = true;
            case 'uniform_log'
                nonlinearity.is_uniform_log = true;
            case 'adapted_log'
                nonlinearity.is_adapted_log = true;
            otherwise
                nonlinearity.is_custom = true;
        end
    else
        nonlinearity.is_custom = true;
    end
else
    nonlinearity = struct('is_modulus',true);
end
if isfield(nonlinearity,'is_modulus') && nonlinearity.is_modulus
    found = true;
else
    nonlinearity.is_modulus = false;
end
if isfield(nonlinearity,'is_uniform_log') && nonlinearity.is_uniform_log
    if found
        error('Parsing error');
    end
    found = true;
    nonlinearity.denominator = default(nonlinearity,'denominator',1.0);
else
    nonlinearity.is_uniform_log = false;
end
if isfield(nonlinearity,'is_adapted_log') && nonlinearity.is_adapted_log
    if found
        error('Parsing error');
    end
    found = true;
    nonlinearity.denominator = default(nonlinearity,'denominator',1.0);
end
if isfield(nonlinearity,'is_custom') && nonlinearity.is_custom
    if found
        error('Parsing error');
    end
else
    nonlinearity.is_custom = false;
end
end

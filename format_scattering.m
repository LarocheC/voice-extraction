function [Smat_1, Smat_2, refs_1, refs_2] = format_scattering(S)
    % Generate references
    spatial_subscripts = 1;
    refs_1 = ...
        generate_refs(S{1+1}.data, spatial_subscripts, S{1+1}.ranges{1+0});
    refs_2 = ...
        generate_refs(S{1+2}.data, spatial_subscripts, S{1+2}.ranges{1+0});

    % Convert scattering nodes to matrices
    nFrames = size(S{1+1}.data, 1);
    nRefs_1 = size(refs_1, 2);
    Smat_1 = zeros(nRefs_1, nFrames);
    for ref1 = 1:nRefs_1
        Smat_1(ref1, :) = subsref(S{1+1}.data, refs_1(ref1));
    end
    nRefs_2 = size(refs_2, 2);
    Smat_2 = zeros(nRefs_2, nFrames);
    for ref2 = 1:nRefs_2
        Smat_2(ref2, :) = subsref(S{1+2}.data, refs_2(:, ref2));
    end

    % Ensure non-negativity
    Smat_1(Smat_1<0)=0;
    Smat_2(Smat_2<0)=0;
end
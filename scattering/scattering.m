function [S, U, Y] = scattering(x, archs, S_log2_oversampling , nmf_log2_oversampling)
    % Compute scattering transform
    [S, U, Y] = sc_propagate(x, archs);

    %% Downsample manually
    phi_oversampling = 2^S_log2_oversampling;
    nmf_oversampling = pow2(nmf_log2_oversampling);
    downsampling = phi_oversampling / nmf_oversampling;
    S{1+1}.data = downsampling * S{1+1}.data(1:downsampling:end, :);
    S{1+1}.ranges{1+0}(2,1) = S{1+1}.ranges{1+0}(2,1) * downsampling;
    archs{1}.banks{1}.behavior.S.log2_oversampling = nmf_log2_oversampling;
    archs{1}.invariants{1}.behavior.S.log2_oversampling = nmf_log2_oversampling;
    for lambda2 = 1:length(S{1+2}.data)
        S{1+2}.data{lambda2} = downsampling * ...
            S{1+2}.data{lambda2}(1:downsampling:end, :);
        S{1+2}.ranges{1+0}{lambda2}(2,1) = ...
            S{1+2}.ranges{1+0}{lambda2}(2,1) * downsampling;
    end
end
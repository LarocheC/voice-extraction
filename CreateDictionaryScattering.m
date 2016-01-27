function W_perc = CreateDictionaryScattering(x, archs, phi_log2_oversampling, nmf_log2_oversampling,option)


    [S, U, Y] = scattering(x, archs, phi_log2_oversampling, nmf_log2_oversampling);
    [Smat_1, Smat_2] = format_scattering(S);
    
    if strcmp(option,'mean')
        W_perc = mean(cat(1, Smat_1, Smat_2),2);
    else
        W_perc = cat(1, Smat_1, Smat_2);
    end


end
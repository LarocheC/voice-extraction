function x_source = mask_scattering( ... 
    S1_source, S2_source, U_mixture, Y_mixture,archs)
Y2_source = dS_backto_dY(S1_source, archs{2});
Y2_source = copy_metadata(Y_mixture{1+1}{1}, Y2_source);
U1fromS1_source = dY_backto_dU(Y2_source);
U1fromS1_source = perform_ft(U1fromS1_source, archs{1}.banks{1}.behavior.key);

Y3_source = dS_backto_dY(S2_source, archs{3});
Y3_source = copy_metadata(Y_mixture{1+2}{1}, Y3_source);
U2_source = dY_backto_dU(Y3_source);

Y2_source = Y_mixture{2}{end};
for lambda2 = 1:length(U2_source.data)
    Y2_source.data{lambda2} = U2_source.data{lambda2} .* ...
        Y_mixture{2}{end}.data{lambda2} ./ U_mixture{1+2}.data{lambda2};
end

Y1_source = dual_scatter_dY(Y2_source, archs{2}.banks{1}, U1fromS1_source);
U1_source = dY_backto_dU(Y1_source);
for lambda1 = 1:length(U1fromS1_source.data)
    Y1_source.data{lambda1} = U1fromS1_source.data{lambda1} + ...
        U1_source.data{lambda1} .* ...
        Y_mixture{1}{end}.data{lambda1} ./ U_mixture{1+1}.data{lambda1};
end

Y0_source = Y_mixture{1+0}{1};
Y0_source.data = complex(zeros(size(Y0_source.data)));
Y0_source.data_ft = complex(zeros(size(Y0_source.data_ft)));
Y0_source = dual_scatter_dY(Y1_source, archs{1}.banks{1}, Y0_source);
x_source = real(Y0_source.data);
end


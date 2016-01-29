function [ W, W2, H2, e] = SPNMF( X , wini , W2, H2 , niter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



%% optimisation de W

W = wini ;
e = zeros(1,niter);
[F,N] = size(X);
XX = X*X';

for iter=1:niter % Main Loop

    W = W .*( (2*XX*W )./ (2*(X*(H2'*W2') + (W2*H2)*X')*W + XX*W*(W'*W) + W*(W'*XX)*W + eps));
    W = W./norm(W);
    
    W2 = W2 .* (X*H2' ) ./ (2*W*(W'*X)*H2' + (W2*H2)*H2' + eps);    
    H2 = H2 .* (W2'*X ) ./ (2*W2'*(W*(W'*X)) + W2'*(W2*H2) +eps) ;  

    % Norm-2 normalization
    scale = sqrt(sum(W2.^2,1)); 
    W2 = W2 .* repmat((scale+eps).^-1,F,1);
    H2 = H2 .* repmat(scale',1,N); 
    
    
      
    e(iter) = norm(X - W*W'*(X) - W2*H2) ;
    
    if mod(iter,100)==0
        display(strcat('Niter',num2str(iter)));
    end
    


end
    
end
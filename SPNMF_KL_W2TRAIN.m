function [ W, H2, e] = SPNMF_KL_W2TRAIN( X , W , W2, H2 , max_iter)
%
% compute SPNMF based on I-divergence (non-nomarlized KL-divergence)
% input:
%   X          nonnegative data input (m times n)
%   max_iter   maximum number of iterations

% output:
%   W          the harmonic dictionary matrix
%   H2         the percussive activation matrix

e=0;
Xsum = sum(X,2);


for iter=1:max_iter

    Z = X ./ (W*(W'*X)+ W2*H2 + eps);
    W = W .* sqrt((Z*(X'*W) +X*(Z'*W)) ...
        ./ (bsxfun(@plus, Xsum'*W, bsxfun(@times, Xsum, sum(W)))+ eps));
    
    H2 = H2 .* ((W2'*X) ./ (W2'*(W*W'*X + W2*H2+eps) +eps)) ;
    
    %e(iter) = idiv(X , W*W'*X + W2*H2);

    if mod(iter,100)==0
        display(strcat('Niter',num2str(iter)));
    end
        
end

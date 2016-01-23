function [ W, W2, H2, e] = SPNMF_KL( X , W , W2, H2 , max_iter)
%
% compute SPNMF based on I-divergence (non-nomarlized KL-divergence)
% input:
%   X          nonnegative data input (m times n, non-sparse)
%   r          number of PNMF components
%   max_iter   maximum number of iterations 

% output:
%   W          the harmonic dictionary matrix
%   W2         the percussice dictionary matrix
%   H2         the percussive activation matrix



e=0;
Xsum = sum(X,2);


[F,T] = size(X);

for iter=1:max_iter
    
    W2 = W2 .* ((X*H2') ./ ((W*W'*X + W2*H2+eps)*H2' +eps)) ;
    H2 = H2 .* ((W2'*X) ./ (W2'*(W*W'*X + W2*H2+eps) +eps)) ;
    
        % Norm-2 normalization
    scale = sqrt(sum(W2.^2,1)); 
    W2 = W2 .* repmat((scale+eps).^-1,F,1);
    H2 = H2 .* repmat(scale',1,T);    
    

    Z = X ./ (W*(W'*X)+ W2*H2 + eps);
    W = W .* sqrt((Z*(X'*W) +X*(Z'*W)) ...
        ./ (bsxfun(@plus, Xsum'*W, bsxfun(@times, Xsum, sum(W)))+ eps));
    
    e(iter) = idiv(X , W*W'*X + W2*H2);
    
    if mod(iter,100)==0
        display(strcat('Niter',num2str(iter)));
    end
        
    
    
end

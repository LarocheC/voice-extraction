function [ W, W2, H2, e] = SPNMF_IS( X , W , W2, H2 , max_iter)
%
% compute SPNMF based on IS-divergence
% input:
%   X          nonnegative data input (m times n, non-sparse)
%   r          number of PNMF components
%   max_iter   maximum number of iterations 

% output:
%   W          the harmonic dictionary matrix
%   W2         the percussice dictionary matrix
%   H2         the percussive activation matrix


e=0;


[F,T] = size(X);
V_ap = W2*H2 + W*W'*X;


for iter=1:max_iter

    
     % Update W
    W2 = W2 .* ((X.*(V_ap+eps).^-2)*H2')./((V_ap+eps).^-1*H2'+eps);
    V_ap = W2*H2 + W*W'*X;
        
    % Update H
    H2 = H2 .* (W2'*(X.*(V_ap+eps).^-2))./(W2'*(V_ap+eps).^-1+eps);
    V_ap = W2*H2 + W*W'*X;

    % Norm-2 normalization
    scale = sqrt(sum(W2.^2,1)); 
    W2 = W2 .* repmat((scale+eps).^-1,F,1);
    H2 = H2 .* repmat(scale',1,T);    


    Z = X ./ (W*(W'*X)+ W2*H2 + eps).^2;
    Z1 = ones(size(X)) ./ (W*(W'*X)+ W2*H2 + eps);
    
    W = W .* ((Z*(X'*W) +X*(Z'*W)) ...
        ./ ((Z1*(X'*W) +X*(Z1'*W))+ eps));
    V_ap = W2*H2 + W*W'*X;

    
% Compute cost value
    e(iter) = sum(X(:)./(V_ap(:)+eps) - log(X(:)./(V_ap(:)+eps)+eps)) - F*T;
    
    if mod(iter,100)==0
        display(strcat('Niter',num2str(iter)));
    end
        
    
    
end

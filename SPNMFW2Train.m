function [ W , H2, e] = SPNMFW2Train( v , wini , w2, hini , niter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



%% optimisation de W

W = wini ;
H2 = hini ;
e = zeros(1,niter);
[F,N] = size(v);
vv = v*v';


for iter=1:niter % Main Loop    
    
    W = W .*  ( (2*vv*W ) ./ (2*(v*(H2'*w2') + (w2*H2)*v')*W + vv*W*(W'*W) + W*(W'*vv)*W  + eps)) ;  
    W = W./norm(W);  
    
	H2 = H2 .* (w2'*v ) ./ (2*w2'*(W*(W'*v)) + w2'*(w2*H2) +eps) ; 


    e(iter) = norm(v - W*W'*(v) - w2*H2,2) ;
    
    if mod(iter,100)==0
        display(strcat('Niter',num2str(iter)));
    end
    


end
    
end


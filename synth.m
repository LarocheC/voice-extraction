function [ s , fp ] = synth( f , t , nbh)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
s = zeros(1,length(t));
i=1;
fp = zeros(1, 15);

for i=1:nbh
    
s = s + 1/(i)^2 * sin (2* pi* i* f* t );

fp(1,i) = i*f; 

end
end


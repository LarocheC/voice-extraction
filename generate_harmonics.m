function [s, fp] = generate_harmonics(f, t, nbh)

s = zeros(1,length(t));
fp = zeros(1, nbh);

for i=1:nbh
    s = s + 1/(i)^2 * sin (2* pi* i* f* t );
    fp(1,i) = i*f;
end
end


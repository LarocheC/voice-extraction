function g = gamma_window(N, alpha, order)
%%
if nargin < 3
    order = 3;
end
if nargin < 2
    alpha = 10.0;
end
t = linspace(0, 1, N);
g = t.^(order-1) .* exp(- alpha * t);
g = g ./ max(g);
end
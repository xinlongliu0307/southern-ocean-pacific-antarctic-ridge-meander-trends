function [h, p, S] = Mann_Kendall(x)
x = x(:); n = length(x);
S = 0;
for k = 1:n-1
    S = S + sum(sign(x((k+1):n) - x(k)));
end
varS = n*(n-1)*(2*n+5)/18;
if S > 0
    Z = (S - 1)/sqrt(varS);
elseif S < 0
    Z = (S + 1)/sqrt(varS);
else
    Z = 0;
end
p = 2 * (1 - normcdf(abs(Z), 0, 1));
h = p < 0.05;
end

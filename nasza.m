function [res] = nasza(n, a, x)
res = 0;
for j=0:n
    res = res + a(j+1)*cos(j*x);
end
end

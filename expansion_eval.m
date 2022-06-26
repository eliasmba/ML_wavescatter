function fx = expansion_eval(u, n, N_multi, centers, x)
%evaluate function with coefficients (u)_m,l in spherical harmonics at x
assert(length(u) == n*N_multi);
assert(length(centers) == n);

Y = zeros(n*N_multi);
for i = 1:n
    [ph, th, ~] = cartsph(x - centers(i));
    
    for l = 1:N_multi+1
        Y(i*(N_multi+1) + l) = harmonicY(0, l, th, ph);
    end
end
fx = u'*Y;
end
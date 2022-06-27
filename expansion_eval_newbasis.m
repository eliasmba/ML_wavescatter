function fx = expansion_eval_newbasis(u, n, N_multi, centers, x, k)
%evaluate function with coefficients (u)_m,l
%in spherical harmonics and bessel basis Yml*H_l(k*.)
%at  point x
assert(length(u) == n*N_multi);
assert(length(centers) == n);

Y = zeros(n*N_multi);
for i = 1:n
    [ph, th, r] = cartsph(x - centers(i));
    
    for l = 1:N_multi+1
        Y(i*(N_multi+1) + l) = 1/sqrt(r)*besselh(l-1+1/2, r* k)*harmonicY(0, l, th, ph);
    end
end
fx = u'*Y;
end
function fx = expansion_eval_newbasis_syms(u, n, N_multi, centers, x, k)
%evaluate function with coefficients (u)_m,l
%in spherical harmonics and bessel basis Yml*H_l(k*.)
%at  point x
assert(length(u) == n*(N_multi + 1));
assert(size(centers,2) == n);
syms Y [n*(N_multi+1) n*(N_multi+1)];
disp(size(centers,1));
for i = 1:n
    disp(x - centers(:,i));
    [ph, th, r] = cart2sph(x(1) - centers(1,i),x(2) - centers(2,i), x(3) - centers(3,i));
    
    for l = 1:N_multi+1
        Y(i*(N_multi+1) + l) = 1/sqrt(r)*besselh(l-1+1/2, r* k)*harmonicY(l, 0, th, ph);
    end
end
fx = u'*Y;
end
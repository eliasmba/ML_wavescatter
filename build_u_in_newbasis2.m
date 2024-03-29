function [u] = build_u_in_newbasis2(sources, intensities,N_multi, centers,k)
%u_in is approixmated by minimizing the L2 norm.
%Output:

n= length(centers);

%syms x y c t v[1 n*(N_multi+1)]

norm_y = @(y) sqrt(y(1,:).^2+y(2,:).^2 + y(3,:).^2);
u_in = @(x) exp(norm_y(repmat(x,1,size(sources, 2)) - sources.'))*intensities;


Y = zeros(n*N_multi);

for i = 1:n
    ph = acos(c(3)/norm_y(c));
    th = atan(c(2)/c(1));
    r = norm_y(c - centers(i));
    
    for l = 1:N_multi+1
        Y(i*(N_multi+1) + l) = 1/sqrt(r(t))*besselh(l-1+1/2, r(t)* k)*harmonicY(0, l, th(t), ph(t));

    end
end
L2_funct = 1;
func = @(v) v*Y;
L2_funct = integrate(L2_funct-func, 0, 1);
u = fminsearch(L2_funct, ones(1, n));
end
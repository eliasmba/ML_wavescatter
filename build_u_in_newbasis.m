function [u] = build_u_in_newbasis(sources, intensities,N_multi, centers,k)
%u_in is approixmated by minimizing the L2 norm.
%Output:

n= length(centers);

syms x y c t v[1 n*N_multi]

norm_y = sqrt(y(1,:).^2+y(2,:).^2 + y(3,:).^2);
u_in = exp(norm_y(repmat(x,1,size(sources, 2)) - sources.'))*intensities;


Y = zeros(n*N_multi);

for i = 1:n
    [ph, th, r] = [acos(c(3)/norm_y(c)), atan(c(2)/c(1)), norm_y(c - centers(i))];
    
    for l = 1:N_multi+1
        Y(i*(N_multi+1) + l) = 1/sqrt(r(t))*besselh(l-1+1/2, r(t)* k)*harmonicY(0, l, th(t), ph(t));

    end
end
L2_funct = integrate(L2_func-v*Y,)
fun
end
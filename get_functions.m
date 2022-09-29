function [M, Sigmax, Sigmaz] = get_functions(sources, receivers, centers, intensities, R, N_multi, kappa0, rho0, kappa_b, rho_b)
%OUTPUT: vector of coefficients of a_n of lenght Nres, vector of eigenmodes
%evaluated at the centers c_r
n_r = size(receivers, 2);
Nres = size(centers, 2);
n_s = size(sources, 2);
omega = 0.01;
k = omega*sqrt(rho_b/kappa_b);

[res, eigenmodes, V] = resonances(centers(1,:), centers(2,:), centers(3,:), R, N_multi, kappa0, rho0, kappa_b, rho_b);

%% Modal Decomposition
om = omega*ones(Nres);
Omega = diag(om - res);
lhs = V*Omega;

vols = 4/3*pi*R.^3;

%u_in = build_u_in_newbasis2(sources, intensities, N_multi, centers,intensities);
u_inc = @(x) inc(x, sources,k);
u_in = make_f_newbasis(u_inc, Nres,centers, N_multi, kappa_b);
Su = makeSmat_newbasis(R, centers, 0, N_multi)\u_in;
S = zeros(Nres, 1);
for  j = 1:Nres
    S(j) = sh_integrate(Su,N_multi, j);
end %Single layer potential inverse integrals

rhs = kappa_b/rho0*S./vols;

a = linsolve(lhs, rhs); %decomposition coefficients

SkD = makeSmat_newbasis(R, centers, k, N_multi);
u_in_pro =  SkD*(SkD\u_in); %solution for lfr system in fourier basis

syms y [3 Nres] x [3 n_s] z [3 n_r]
absDy =  (y(1,:).'- y(1,:)).^2+ (y(2,:).'-y(2,:)).^2 + (y(3,:).'-y(3,:)).^2;
absDx =  (y(1,:).'- x(1,:)).^2+ (y(2,:).'-x(2,:)).^2 + (y(3,:).'-x(3,:)).^2;
absDz =  (y(1,:).'- z(1,:)).^2+ (y(2,:).'-z(2,:)).^2 + (y(3,:).'-z(3,:)).^2;

M = zeros(Nres);
Sigmax = zeros(Nres, n_r);
Sigmaz = zeros(Nres, n_s);

for i=1:Nres
    M(i,i) = 1;
    for j=1:Nres
        if i~=j
            M(i,j) = a(i,j)*expansion_eval_newbasis(eigenmodes(i,:), Nres, N_multi, centers,absDy(i,j) , k);
        end
    end

    for j=1:n_r
        Sigmax(i,j) = expansion_eval_newbasis(u_in_pro, Nres, N_multi, centers, absDx(i,j), k);
    end

    for j=1:n_r
        Sigmaz(i,j) = expansion_eval_newbasis(u_in_pro, Nres, N_multi, centers, absDz(i,j), k);
    end
end

end

function [u] = inc(x, sources,k)
    norm_y = @(y) sqrt(y(1,:).^2+y(2,:).^2 + y(3,:).^2);
    u = sum(exp(k*norm_y(repmat(x(1),1,size(sources, 2)) - sources(1,:).')),2);
end
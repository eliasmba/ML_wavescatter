function [z] = LFRsolver(sources, receivers ,intensities, centers, R, k, k0)
%resonances and eigenmodes of the system
n_r = length(receivers);
Nres = length(centers);
[res, eigenmodes, V] = resonances(centers(1), centers(2), centers(3), R, N_multi, kappa0, rho0, kappa_b, rho_b);
%% Modal Decomposition
om = omega*ones(Nres);
Omega = diag(om - res);
lhs = V*Omega;

vols = 4/3*pi*R.^3;

u_in = built_u_in_newbasis(sources, intensities);
Su = makeSmat_newbasis(R, centers, 0, N_multi)\u_in;
S = zeros(Nres, 1);
for  j = 1:Nres
    S(j) = sh_integrate(Su, Nres, N_multi, j, centers(j), R);
end %Single layer potential inverse integrals

rhs = kappa_b/rho0*S./vols;

a = linsolve(lhs, rhs); %decomposition coefficients

SkD = makeSmat_newbasis(R, centers, k, N_multi);
u_in_pro =  SkD*(SkD\u_in); %solution for lfr system in fourier basis
                                         %todo: u_in
z = zeros(n_r);                                  
for i = 1:n_r
    z(i) = expansion_eval_newbasis(u_in_pro, Nres, N_multi, centers, receivers(i), k); %u_in
    for k = 1:Nres
        z(i) = z(i) + a(k)*expansion_eval(eigenmodes(k,:), Nres, N_multi, centers, receivers(i), k);
    end
end
end

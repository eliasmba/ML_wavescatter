function [u] = LFRsolver()
%resonances and eigenmodes of the system
[res, eigenmodes, V] = resonances(cx, cy, cz, R, N_multi, kappa0, rho0, kappa_b, rho_b);
%% Modal Decomposition
om = omega*ones(Nres);
Omega = diag(om - res);
lhs = V*Omega;

vols = 4/3*pi*R.^3;
S = %Single layer potential inverse integrals ???
rhs = kappa_b/rho0*S./vols;

an = linsolve(lhs, rhs); %decomposition coefficients

u = sum(eigenmodes.*an); %solution for lfr system in fourier basis
                                         %todo: u_in


end
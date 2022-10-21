function [resonances, eigenmodes, V] = resonances(cx, cy, cz, R, N_multi, kappa0, rho0, kappa_b, rho_b)
%returns eigenmodes and resonances in fourier basis for circular low frequency resonators
%V passed for modal decomposition
%% THM 2.7
delta = rho0/rho_b;
Nres = length(R);
C = capacitance(cx, cy, cz ,R,rho0,rho_b,kappa0,kappa_b,delta);
omega = 0.01;
k = omega*sqrt(rho_b/kappa_b);

[V, resonances] = eig(C);
resonances = sqrt(resonances);
centers = [cx;cy;cz];
S0 = MakeSmat_newbasis(R(1), centers, 0.00001, 0); %here we set N_multi = 0, as the characteristic function is constant on the resonators(we therefore only need m=l=0)

psi = sqrt(4)*S0^-1; %we have psi1,..., psiN in the columns
%% Cor 2.9
Sk = MakeSmat_newbasis(R(1), centers, k, 0); %Single layer potential for spherical resonators in spherical harmonics

Sk_psi = (Sk*psi).'; %= SkD(psi1), ...,SkD(psiN) in the rows
eigenmodes = V*Sk_psi';% = Sk_psi*v1, ..., Sk_psi*vn in the rows
                                   %elt of lR^(Nres x Nres(N_multi+1))
                                   %low frequency resonance modes outside of resonators
end

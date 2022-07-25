function [resonances, eigenmodes, V] = resonances(cx, cy, cz, R, N_multi, kappa0, rho0, kappa_b, rho_b)
%returns eigenmodes and resonances in fourier basis for circular low frequency resonators
%V passed for modal decomposition
%% THM 2.7
Nres = length(R);
C = capacitance(cx, cy, cz ,R,rho0,rho_b,kappa0,kappa_b,delta);

[V, resonances] = eig(C);
resonances = sqrt(resonances);

S0 = MakeSmat_newbasis(R(1), centers, 0, N_multi);
psi = zeros(Nres, Nres*(N_multi+1));
for j = 1:Nres
    chi_j = make_char_fct_newbasis(centers(j), R(1));
    psi(j,:) = S0\chi_j;
end
%% Cor 2.9
Sk = MakeSmat_newbasis(R(1), centers, k, N_multi); %Single layer potential for spherical resonators in spherical harmonics

Sk_psi = psi*Sk'; %= SkD(psi1), ...,SkD(psiN) in the rows
eigenmodes = V*Sk_psi';% = Sk_psi*v1, ..., Sk_psi*vn in the rows
                                   %elt of lR^(Nres x Nres(N_multi+1))
                                   %low frequency resonance modes outside of resonators
end

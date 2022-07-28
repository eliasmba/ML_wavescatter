function [resonances, eigenmodes, V] = resonances(cx, cy, cz, R, N_multi, kappa0, rho0, kappa_b, rho_b)
%returns eigenmodes and resonances in fourier basis for circular low frequency resonators
%V passed for modal decomposition
%% THM 2.7
Nres = length(R);
C = capacitance(cx, cy, cz ,R,rho0,rho_b,kappa0,kappa_b,delta);

[V, resonances] = eig(C);
resonances = sqrt(resonances);

S0 = MakeSmat_newbasis(R, centers, 0, N_multi);
psi = zeros(Nres, Nres*(N_multi+1));
for j = 1:Nres
    char_j = @(x) (char_f_sphere(x, centers(j), R));
    chi_j = make_f_newbasis(char_j, N_multi);
    psi(j,:) = S0\chi_j;
end
%% Cor 2.9
Sk = MakeSmat(R, centers, k, N_multi); %Single layer potential for spherical resonators in spherical harmonics

Sk_psi = psi*Sk'; %= SkD(psi1), ...,SkD(psiN) in the rows
eigenmodes = V*Sk_psi';% = Sk_psi*v1, ..., Sk_psi*vn in the rows
                                   %elt of lR^(Nres x Nres(N_multi+1))
                                   %low frequency resonance modes outside of resonators
end

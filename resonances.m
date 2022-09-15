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
S0 = MakeSmat_newbasis(R(1), centers, 0, N_multi);
psi = zeros(Nres, Nres*(N_multi+1));
for j = 1:Nres
    char_j = @(x) (char_f_sphere(x, centers(:,j), R(1)));
    chi_j = make_f_newbasis(char_j, Nres, centers, N_multi, k);
    disp(size(chi_j));
    chi_j = chi_j.';
    psi(j,:) = S0\chi_j;
end
%% Cor 2.9
Sk = MakeSmat_newbasis(R(1), centers, k, N_multi); %Single layer potential for spherical resonators in spherical harmonics

Sk_psi = psi*Sk'; %= SkD(psi1), ...,SkD(psiN) in the rows
eigenmodes = V*Sk_psi';% = Sk_psi*v1, ..., Sk_psi*vn in the rows
                                   %elt of lR^(Nres x Nres(N_multi+1))
                                   %low frequency resonance modes outside of resonators
end

function [ind] = char_f_sphere(x, c, r)
    if abs(norm(x - c) - r) < 0.0001
        ind = 1;
    else
        ind = 0;
   end
end

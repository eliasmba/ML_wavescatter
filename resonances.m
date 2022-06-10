function [resonances, eigenmodes, V] = resonances(cx, cy, cz, R, N_multi, kappa0, rho0, kappa_b, rho_b)
%returns eigenmodes and resonances in fourier basis for circular low frequency resonators
%V passed for modal decomposition
%% THM 2.7
Nres = length(R);
C = capacitance(cx, cy, cz ,R,rho0,rho_b,kappa0,kappa_b,delta);

[V, resonances] = eig(C);
resonances = sqrt(resonances);

psi = stuffidk;
%% Cor 2.9
A = MakeA(R,omega,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy);
SkD = A(1:Nres*(N_multi+1), 1:Nres*(N_multi+1)); %Single layer potential for circular resonators in fourier basis

SkDx = psi*SkD'; %= SkD(psi1), ...,SkD(psiN) in the rows
eigenmodes = V*SkDx';% = SkDx*v1, ..., SkDx*vn in the rows
                                   %elt of lR^(Nres x Nres(N_multi+1))
                                   %low frequency resonance modes outside of resonators
end
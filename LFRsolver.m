function [z] = LFRsolver(sources, receivers ,intensities, centers, R, k, k0)
%resonances and eigenmodes of the system
n_r = size(receivers, 2);
Nres = size(centers, 2);
%n_s = size(sources, 2);
omega = 0.01;
k = omega*sqrt(rho_b/kappa_b);

N_multi=0;

[res, eigenmodes, V] = resonances(centers(1,:), centers(2,:), centers(3,:), R, N_multi, kappa0, rho0, kappa_b, rho_b);


%% Modal Decomposition
om = omega*ones(Nres,1);
Omega = om - res;

lhs = V*Omega;



vols = 4/3*pi*R.^3;

%u_in = build_u_in_newbasis2(sources, intensities, N_multi, centers,intensities);
u_in = inc(centers, sources,3).';

S = MakeSmat_newbasis(R(1), centers, 0.0001, N_multi)\u_in*sqrt(4);
%Single layer potential inverse integrals

rhs = kappa_b/rho0*S./vols;


a = linsolve(lhs, rhs); %decomposition coefficients


SkD = MakeSmat_newbasis(R(1), centers, k, N_multi);
u_in_pro =  SkD*(SkD\u_in); %solution for lfr system in fourier basis
z = zeros(n_r);                                  
for i = 1:n_r
    z(i) = inc(receivers(i), soruces, k); %u_in
    for k = 1:Nres
        z(i) = z(i) + a(k)*eigenmodes(k,:);
    end
end
end

function[u] = inc(x, sources, k)
norm_y = @(y, sources) sqrt((y(1,:)).^2+(y(2,:)).^2 + (y(3,:)).^2);
u = exp(k*norm_y(x,sources));
end

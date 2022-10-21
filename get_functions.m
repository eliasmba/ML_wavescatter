function [M, Sigmax, Sigmaz] = get_functions(sources, receivers, centers, intensities, R, N_multi, kappa0, rho0, kappa_b, rho_b)
%OUTPUT: vector of coefficients of a_n of lenght Nres, vector of eigenmodes
%evaluated at the centers c_r
n_r = size(receivers, 2);
Nres = size(centers, 2);
n_s = size(sources, 2);
omega = 0.01;
k = omega*sqrt(rho_b/kappa_b);

N_multi=0;

[res, eigenmodes, V] = resonances(centers(1,:), centers(2,:), centers(3,:), R, N_multi, kappa0, rho0, kappa_b, rho_b);


%% Modal Decomposition
om = omega*ones(Nres);
Omega = diag(om - res);
lhs = V*Omega;

vols = 4/3*pi*R.^3;

%u_in = build_u_in_newbasis2(sources, intensities, N_multi, centers,intensities);
u_in = inc(centers, sources,k);
disp(size(u_in));
S = MakeSmat_newbasis(R(1), centers, 0.0001, N_multi)\u_in;
%Single layer potential inverse integrals

rhs = kappa_b/rho0*S./vols;

a = linsolve(lhs, rhs); %decomposition coefficients

disp(size(a));

SkD = MakeSmat_newbasis(R(1), centers, k, N_multi);
u_in_pro =  SkD*(SkD\u_in); %solution for lfr system in fourier basis

%syms y [3 Nres] x [3 n_s] z [3 n_r]
x = sources;
y = receivers;
z = centers;

absDy =  (y(1,:).'- y(1,:)).^2+ (y(2,:).'-y(2,:)).^2 + (y(3,:).'-y(3,:)).^2;
absDx =  (y(1,:).'- x(1,:)).^2+ (y(2,:).'-x(2,:)).^2 + (y(3,:).'-x(3,:)).^2;
absDz =  (y(1,:).'- z(1,:)).^2+ (y(2,:).'-z(2,:)).^2 + (y(3,:).'-z(3,:)).^2;

M = zeros(Nres);
Sigmax = zeros(Nres, n_r);
Sigmaz = zeros(Nres, n_s);
disp(Nres);
for i=1:Nres
    M(i,i) = 1;
    for j=1:Nres
        disp(j);
        if i~=j
            M(i,j) = a(i)*eigenmodes(i,:)*(y(i,:)-y(j,:)).'; %this is u_n
        end
    end

    for j=1:n_r
        Sigmax(i,j) = exp(k*absDx(i,j));
    end

    for j=1:n_r
        Sigmaz(i,j) = exp(k*absDz(i,j));
    end
end

end

%function [u] = inc(x, sources,k)
 %   norm_y = @(y) sqrt(y(1,:).^2+y(2,:).^2 + y(3,:).^2);
  %  u = sum(exp(k*norm_y(repmat(x(1),1,size(sources, 2)) - sources(1,:).')),2);
%end
function[u] = inc(x, sources, k)
norm_y = @(y, sources) sqrt((y(1,:).'-sources(1,:)).^2+(y(2,:).'-sources(2,:)).^2 + (y(3,:).'-sources(3,:)).^2);
u = exp(k*norm_y(x,sources));
end
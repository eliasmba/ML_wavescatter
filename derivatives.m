function [D1M, D2M, D3M, D1Sigmax, D2Sigmax, D3Sigmax, D1Sigmaz, D2Sigmaz, D3Sigmaz] = derivatives(sources, receivers, centers, intensities, R, N_multi, kappa0, rho0, kappa_b, rho_b)
n_r = size(receivers, 2);
Nres = size(centers, 2);
n_s = size(sources, 2);
omega = 0.01;
k = omega*sqrt(rho_b/kappa_b);

N_multi=0;

[res, ~, V] = resonances(centers(1,:), centers(2,:), centers(3,:), R, N_multi, kappa0, rho0, kappa_b, rho_b);
eig_diff = eigenmodes(centers(1,:), centers(2,:), centers(3,:), R, N_multi, kappa0, rho0, kappa_b, rho_b);


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
y = centers;
z = receivers;

absDy =  (y(1,:).'- y(1,:)).^2+ (y(2,:).'-y(2,:)).^2 + (y(3,:).'-y(3,:)).^2;
absDx =  (y(1,:).'- x(1,:)).^2+ (y(2,:).'-x(2,:)).^2 + (y(3,:).'-x(3,:)).^2;
absDz =  (y(1,:).'- z(1,:)).^2+ (y(2,:).'-z(2,:)).^2 + (y(3,:).'-z(3,:)).^2;

D1M = zeros(Nres);
D2M = zeros(Nres);
D3M = zeros(Nres);
D1Sigmax = zeros(Nres, n_r);
D2Sigmax = zeros(Nres, n_r);
D3Sigmax = zeros(Nres, n_r);
D1Sigmaz = zeros(Nres, n_s);
D2Sigmaz = zeros(Nres, n_s);
D3Sigmaz = zeros(Nres, n_s);

for i=1:Nres
    D1M(i,i) = 0;
    D2M(i,i) = 0;
    for j=1:Nres
        if i~=j
            
            D1M(i,j) = a(i)*sum(eig_diff(i,:));%*absDy(i,j); %this is u_n
            D2M(i,j) = a(i)*sum(eig_diff(i,:));%*absDy(i,j); %this is u_n
            D3M(i,j) = a(i)*sum(eig_diff(i,:));%*absDy(i,j); %this is u_n
        end
    end

    for j=1:n_r
        D1Sigmax(i,j) = exp(k*absDx(i,j)).*(y(1,i)-y(1,j))./absDy(i,j);
        D2Sigmax(i,j) = exp(k*absDx(i,j)).*(y(2,i)-y(2,j))./absDy(i,j);
        D3Sigmax(i,j) = exp(k*absDx(i,j)).*(y(3,i)-y(3,j))./absDy(i,j);
    end

    for j=1:n_r
        D1Sigmaz(i,j) = exp(k*absDz(i,j)).*(y(1,i)-y(1,j))./absDy(i,j);
        D2Sigmaz(i,j) = exp(k*absDz(i,j)).*(y(2,i)-y(2,j))./absDy(i,j);
        D3Sigmaz(i,j) = exp(k*absDz(i,j)).*(y(3,i)-y(3,j))./absDy(i,j);
    end 
end

D1M(isnan(D1M) | isinf(D1M)) = 0;
D2M(isnan(D2M) | isinf(D2M)) = 0;
D3M(isnan(D3M) | isinf(D3M)) = 0;

D1Sigmax(isnan(D1Sigmax) | isinf(D1Sigmax)) =0;
D2Sigmax(isnan(D2Sigmax) | isinf(D2Sigmax)) =0;
D3Sigmax(isnan(D3Sigmax) | isinf(D3Sigmax)) =0;
D1Sigmaz(isnan(D1Sigmaz) | isinf(D1Sigmaz)) =0;
D2Sigmaz(isnan(D2Sigmaz) | isinf(D2Sigmaz)) =0;
D3Sigmaz(isnan(D3Sigmaz) | isinf(D3Sigmaz)) =0;


end

%function [u] = inc(x, sources,k)
 %   norm_y = @(y) sqrt(y(1,:).^2+y(2,:).^2 + y(3,:).^2);
  %  u = sum(exp(k*norm_y(repmat(x(1),1,size(sources, 2)) - sources(1,:).')),2);
%end
function[u] = inc(x, sources, k)
norm_y = @(y, sources) sqrt((y(1,:).'-sources(1,:)).^2+(y(2,:).'-sources(2,:)).^2 + (y(3,:).'-sources(3,:)).^2);
u = exp(k*norm_y(x,sources));
end
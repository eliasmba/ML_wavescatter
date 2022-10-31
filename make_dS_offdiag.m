function dSk = make_dS_offdiag(ci, cj, R, k, N_multi)
r_x = norm(cj-ci);

Jdata = zeros(N_multi+1,1);
deriJdata = zeros(N_multi+1,1);
J1data = zeros(N_multi+1,1);

deriJdata(1) = -sqrt(pi/2/R)*besselj(3/2, k*r_x);

for l = 1:N_multi+1
    Jdata(l) = besselj(l-1+1/2, k*r_x);
    J1data(l) = besselh(l-1+1/2, k*R);
end

if N_multi > 0
    deriJdata(2) = sqrt(pi/2/(k*r_x))*besselj(3/2,k*r_x);

    for n = 2:N_multi+1
        deriJdata(n) = Jdata(n-1) - Jdata(n)*n/(k*r_x);
    end
end

dSk = zeros(N_multi);
const=-1i*k*R;

for l=1:N_multi+1
    for lp=1:N_multi+1

       Jlpk0= deriJdata(lp);
       Jlk0= J1data(l);

       dSk(l,lp)=const*Jlpk0*A_coeff(lp-1,l-1,N_multi,k*r_x)*Jlk0;

    end
end
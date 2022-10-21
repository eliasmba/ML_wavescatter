function Sk0 = makeS_offdiag(ci,cj,R, k, N_multi)
% make the off-diagonal blocks of S(r)
r_x = norm(cj-ci);

Jdata = zeros(N_multi+1,1);
J1data = zeros(N_multi+1,1);

for l = 1:N_multi+1
    Jdata(l) = besselj(l-1+1/2, k*r_x);
    J1data(l) = besselh(l-1+1/2, k*R);
end

Sk0=zeros(N_multi+1);

const=-1i*k*R;

for l=1:N_multi+1
    for lp=1:N_multi+1

       Jlpk0= Jdata(lp);
       Jlk0= J1data(l);

       Sk0(l,lp)=const*Jlpk0*A_coeff(lp-1,l-1,N_multi,k*r_x)*Jlk0;

    end
end
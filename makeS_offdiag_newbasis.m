function Sk0 = makeS_offdiag_newbasis(ci,cj,R, k, N_multi)
% make the off-diagonal blocks of S(r)
r_x = norm(cj-ci);

Jdata = zeros(N_multi+1,1);
J1data = zeros(N_multi+1,1);
Hdata = zeros(N_multi+1, 1);

for l = 1:N_multi+1
    Jdata(l) = besselj(l-1+1/2, k*r_x);
    Hdata(l) = besselh(l-1+1/2, k*r_x);
    J1data(l) = besselh(l-1+1/2, k*R);
end

Sk=zeros(N_multi+1);

const=-1i*pi/2*sqrt(R^3);

for l=1:N_multi+1
    for lp=1:N_multi+1

       Jlpk= Jdata(lp);
       Jlk= J1data(l);
       Hlk = Hdata(l);

       Sk(l,lp)=const*Jlpk*A_coeff(lp-1,l-1,N_multi,k*r_x)*Jlk/Hlk;

    end
end
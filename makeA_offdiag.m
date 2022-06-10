function M = makeA_offdiag(i,j,ci,cj,Ri,Rj,k0,delta,N_multi,SpecialFuncDataLi,SpecialFuncDataLj)
% make the off-diagonal blocks of A


%SpecialFuncDataL=[Jdata_k0R,dJdata_k0R];
Jdata_k0Ri=SpecialFuncDataLi(:,1);
dJdata_k0Ri=SpecialFuncDataLi(:,2);

Jdata_k0Rj=SpecialFuncDataLj(:,1);

Sk0=zeros(N_multi+1);
dSk0=zeros(N_multi+1);
zeromatrix=zeros(N_multi+1);
const=-1i*k0*Rj^2;

r_x = norm(cj-ci);

for l=1:N_multi+1
for lp=1:N_multi+1
   
   %Jnk0= besselj(n,k0*R);
   %Jmk0= besselj(m,k0*R);
   %Hnmk0= besselh(n-m,1,k0*mag_cicj);
   %dJmk0= 1/2*(besselj(m-1, k0*R) - besselj(m+1, k0*R));
   Jlpk0= Jdata_k0Rj(lp);
   Jlk0= Jdata_k0Ri(l);
   dJlk0= dJdata_k0Ri(l);
   
%    if j>=i
%     Hnmk0= Hdata_cicj(i,j,N_multi2+1+n-m);
%    else
%     Hnmk0= Hdata_cicj(j,i,N_multi2+1+n-m);
%    end
   
   Sk0(l,lp)=const*Jlpk0*A_coeff(lp-1,l-1,N_multi,k0*r_x)*Jlk0;
   dSk0(l,lp)=const*k0*Jlpk0*A_coeff(lp-1,l-1,N_multi,k0*r_x)*dJlk0;
   
end
end

M=[zeromatrix, -Sk0; zeromatrix, -delta*dSk0];

end
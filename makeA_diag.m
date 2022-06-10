function M = makeA_diag(R,k0,kb,delta,N_multi,SpecialFuncDataM)
% make the diagonal blocks of A

Jdata_k0R=SpecialFuncDataM(:,1);
Jdata_kbR=SpecialFuncDataM(:,2);
Hdata_k0R=SpecialFuncDataM(:,3);
Hdata_kbR=SpecialFuncDataM(:,4);
dJdata_kbR=SpecialFuncDataM(:,5);
dHdata_k0R=SpecialFuncDataM(:,6);

const = -1i*R^2;

Sk0=zeros(N_multi+1);
Skb=zeros(N_multi+1);
dSk0=zeros(N_multi+1);
dSkb=zeros(N_multi+1);

for n=1:N_multi+1
   
   Jk0 = Jdata_k0R(n);
   Jkb = Jdata_kbR(n);
   Hk0 = Hdata_k0R(n);
   Hkb = Hdata_kbR(n);
   dHk0 = dHdata_k0R(n);
   dJkb = dJdata_kbR(n);
    
   Sk0(n,n) = const*k0*Jk0*Hk0;
   dSk0(n,n) = const*k0^2*Jk0*dHk0;
   Skb(n,n) = const*kb*Jkb*Hkb;
   dSkb(n,n) = const*kb^2*Hkb*dJkb;
   
end
M=[Skb, -Sk0; dSkb, -delta*dSk0];

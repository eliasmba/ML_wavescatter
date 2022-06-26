function Skb = makeS_diag(R, r, k,N_multi)
% make the diagonal blocks of S

%Jdata_kbR = makeBesselJdata(N_multi, k*R) =SpecialFuncDataM(:,2);

Jdata = zeros(N_multi+1,1);
Hdata = zeros(N_multi+1,1);

for l=1:N_multi+1
    Jdata(l)= besselj(l-1+1/2, k*R);
    Hdata(l)= besselh(l-1+1/2, k*r);
end

const = -1i*pi/2*sqrt(R^3/r);
Skb=zeros(N_multi+1);

for l=1:N_multi+1
   Jk = Jdata(l);
   Hk = Hdata(l);
    
   Skb(l,l) = const*Jk*Hk;
   
end


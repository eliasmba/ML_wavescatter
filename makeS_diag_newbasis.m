function Skb = makeS_diag_newbasis(R, k,N_multi)
% make the diagonal blocks of S

const = -1i*pi/2*sqrt(R^3);
Skb=zeros(N_multi+1);

for l=1:N_multi+1
   Jk = besselj(l-1+1/2, k*R);
    
   Skb(l,l) = const*Jk;
   
end


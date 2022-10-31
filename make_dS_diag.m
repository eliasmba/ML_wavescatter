function dSkb = make_dS_diag(R, r, k,N_multi)
% make the diagonal blocks of S

%Jdata_kbR = makeBesselJdata(N_multi, k*R) =SpecialFuncDataM(:,2);

Jdata =  zeros(N_multi+1, 1);
deriJdata = zeros(N_multi+1,1);
Hdata = zeros(N_multi+1,1);
deriJdata(1) = -sqrt(pi/2/R)*besselj(3/2, k*R);

for l=1:N_multi+1
    Jdata(l)= sqrt(pi/2/R)*besselj(l-1+1/2, k*R);
    Hdata(l)= sqrt(pi/2/r)*besselh(l-1+1/2, k*r);
end

if N_multi > 0
    deriJdata(2) = sqrt(pi/2/(k*R))*besselj(3/2,k*R);

    for n = 2:N_multi+1
        deriJdata(n) = Jdata(n-1) - Jdata(n)*n/(k*R);
    end
end

const = -1i*k^2*R^2;
dSkb=zeros(N_multi+1);

for l=1:N_multi+1
   Jk = Jdata(l);
   Hk = Hdata(l);
    
   dSkb(l,l) = const*Jk*Hk*k;
   
end


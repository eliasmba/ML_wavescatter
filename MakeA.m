function A = MakeA(R,omega,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy)

N = length(cx);

k0=omega*sqrt(rho0/kappa0);
kb=omega*sqrt(rho_b/kappa_b);

Jdata_k0 = [];
Jdata_kb = [];
Hdata_k0 = [];
Hdata_kb = [];
dJdata_k0 = [];
dJdata_kb = [];
dHdata_k0 = [];

for i = 1:N
Jdata_k0 = [Jdata_k0, makeBesselJdata(N_multi, k0*R(i)  )];
Jdata_kb = [Jdata_kb, makeBesselJdata(N_multi, kb*R(i)  )];
Hdata_k0 = [Hdata_k0, makeHankel1data(N_multi, k0*R(i)  )];
Hdata_kb = [Hdata_kb, makeHankel1data(N_multi, kb*R(i)  )];
dJdata_k0 = [dJdata_k0, makeDeriBesselJdata(N_multi,k0*R(i),Jdata_k0(:,i))];
dJdata_kb = [dJdata_kb, makeDeriBesselJdata(N_multi,kb*R(i),Jdata_kb(:,i))];
dHdata_k0 = [dHdata_k0, makeDeriHankel1data(N_multi,k0*R(i),Hdata_k0(:,i))];
end


%%


SpecialFuncDataM = [];
SpecialFuncDataL = [];
for i = 1:N
SpecialFuncDataM = [SpecialFuncDataM, Jdata_k0(:,i),Jdata_kb(:,i),Hdata_k0(:,i),Hdata_kb(:,i),dJdata_kb(:,i),dHdata_k0(:,i)];
SpecialFuncDataL =[SpecialFuncDataL, Jdata_k0(:,i),dJdata_k0(:,i)];
end


N_oneblock = (N_multi+1)*2; 

A = zeros(N_oneblock*N);

for i=1:N
    for j=1:N
        if i==j
            A((i-1)*N_oneblock+1:i*N_oneblock,(i-1)*N_oneblock+1:i*N_oneblock)= makeA_diag(R(i),k0,kb,delta,N_multi,SpecialFuncDataM(:,6*(i-1)+1:6*i));
        else
            A((i-1)*N_oneblock+1:i*N_oneblock,(j-1)*N_oneblock+1:j*N_oneblock)= makeA_offdiag(i,j,[cx(i),cy(i)],[cx(j),cy(j)],R(j),R(i),k0,delta,N_multi,SpecialFuncDataL(:,2*(i-1)+1:2*i),SpecialFuncDataL(:,2*(j-1)+1:2*j));
        end
    end
end

end
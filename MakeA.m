function A = MakeA(R,omega,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy,cz)

N = length(cx);
R = R(1);

k0=omega*sqrt(rho0/kappa0);
kb=omega*sqrt(rho_b/kappa_b);

Jk0 = sqrt(pi/2/R/k0)*besselj(1/2, k0*R);
Jkb = sqrt(pi/2/R/kb)*besselj(1/2,kb*R);
Hk0 = sqrt(pi/2/R/k0)*besselh(1/2,1,k0*R);
Hkb = sqrt(pi/2/R/kb)*besselh(1/2,1,kb*R);
dJk0 = -sqrt(pi/2/R/k0)*besselj(3/2,k0*R);
dJkb = -sqrt(pi/2/R/kb)*besselh(3/2,kb*R);
dHk0 = -sqrt(pi/2/R/k0)*besselh(3/2,1,k0*R);

const = -1i*R^2;

Sk0 = const*k0*Jk0*Hk0;
dSk0 = const*k0^2*Jk0*dHk0;
Skb = const*kb*Jkb*Hkb;
dSkb = const*kb^2*Hkb*dJkb;

M=[Skb, -Sk0; dSkb, -delta*dSk0];


const2=-1i*k0*R^2;
r = sqrt((cx.'-cx)^2 + (cy.'-cy)^2 + (cz.'-cz)^2); 
coeff_A = @(y)sqrt(1/4/pi)*C_coeff(0,0,0,0,0,0)*sqrt(pi/2/y)*besselh(1/2,1,y);

A = zeros(2*N);

for i=1:N
    for j=1:N
        if i==j
            A(2*(i-1)+1:2*i,2*(i-1)+1:2*i)= M;
        else
            A(2*(i-1)+1:2*i,2*(j-1)+1:2*j)= [0, -const2*Jk0*coeff_A(k0*r(i,j))*Jk0;0,const2*k0*Jk0*coeff_A(k0*r(i,j))*dJk0];
        end
    end
end

end
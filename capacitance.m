function C = capacitance(cx, cy, cz ,R,rho0,rho_b,kappa0,kappa_b,delta)
%TODO: fit for cx, cy, cz & R a constant scalar
N = length(cx);

N_multi = 0;

omega = 0.000001;
A = MakeA(R,omega,rho0,rho_b,kappa0,kappa_b,delta,N_multi,cx,cy);

B = zeros(N);
for i = 1:N
    B(i,i) = A(2*i-1,2*i-1);
end
for i = 1:N
    for j = 1:N
        if i ~= j
            B(i,j) = -A(2*i-1,2*j);
        end
    end
end

B = -B/4/pi/R(1)^2;

C = B^-1;
end

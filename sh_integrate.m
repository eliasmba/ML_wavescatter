function [I] = sh_integrate(u,N_multi, i)

I =u(i*N_multi, i*N_multi)*sqrt(4*pi);
end
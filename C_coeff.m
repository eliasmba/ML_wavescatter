function C = C_coeff(lp, mp, l, m, lambda, mu)

C = 1i^(l-lp+lambda)*(-1)^mp*sqrt(4*pi*(2*l+1)*(2*lp+1)*(2*lambda+1))*...
    Wigner3j([lp l lambda],[0 0 0])*Wigner3j([lp l lambda],[-mp m mu]);

end
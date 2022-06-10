function A = A_coeff(lp,l,N_multi,z)

    Hdata = makeHankel1data(N_multi,z);

    A = 0;
    for n = 0:N_multi
        A = A + sqrt((2*n+1)/4/pi)*C_coeff(lp,0,l,0,n,0)*Hdata(n+1);
    end
    
end

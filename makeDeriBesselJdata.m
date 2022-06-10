function deriJdata = makeDeriBesselJdata(N_multi,z,Jdata)
% Computes the derivative of the Bessel function j_l, evaluated at z, using
% recurrence relations

if N_multi == 0
    deriJdata = -sqrt(pi/2/z)*besselj(3/2,z);
else
    deriJdata=zeros(N_multi+1,1);
    
    deriJdata(1) = -Jdata(2);
    
    for n = 2:N_multi+1
        deriJdata(n) = Jdata(n-1) - Jdata(n)*n/z;
    end

%     n = -N_multi;
%     deriJdata(1) = n/z*Jdata(N_multi+1+n) - Jdata(N_multi+1+n+1)  ;

end

end
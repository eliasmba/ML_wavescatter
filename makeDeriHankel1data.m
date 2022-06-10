function deriHdata = makeDeriHankel1data(N_multi,z,Hdata)
% Computes the derivative of the Hankel function h_l^(1), evaluated at z,
% using recurrence relations

if N_multi == 0
    deriHdata = -sqrt(pi/2/z)*besselh(3/2,1,z);
else
    deriHdata=zeros(N_multi+1,1);
    
    deriHdata(1) = -Hdata(2);
    
    for n = 2:N_multi+1
        deriHdata(n) = Hdata(n-1) - Hdata(n)*n/z;
    end
    
    
%     for n = -(N_multi-1):N_multi
%         deriHdata(N_multi+1+n) = Hdata(N_multi+1+n-1) - Hdata(N_multi+1+n)*(n+1)/z;
%     end

%     n = -N_multi;
%     deriHdata(1) = n/z*Hdata(N_multi+1+n) - Hdata(N_multi+1+n+1)  ;

end

end
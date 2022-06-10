function Jdata = makeBesselJdata(N_multi,z)
% Gives the spherical Bessel function j_l evaluated at z. The first two
% (l=0, l=1) are defined explicitly then a recurrence relation is used to
% generate the remaining l=0,1,...,N_multi

Jdata = zeros(N_multi+1,1);

Jdata(1) = sqrt(pi/2/z)*besselj(1/2,z);
if N_multi > 0
    Jdata(2) = sqrt(pi/2/z)*besselj(3/2,z);

    for n=3:N_multi+1
        Jdata(n)=-Jdata(n-2)+(2*n-3)/z*Jdata(n-1);
    end
end

end


% Jdata = zeros(2*N_multi+1,1);
% 
% Jdata(N_multi+1+0) = sqrt(pi/2/z)*besselj(1/2,z);
% Jdata(N_multi+1+1) = sqrt(pi/2/z)*besselj(3/2,z);
% 
% for n=2:N_multi
%     Jdata(N_multi+1+n)=-Jdata(N_multi+1+n-2)+(2*n-1)/z*Jdata(N_multi+1+n-1);
% end
% 
% for n=-1:-1:-N_multi
%     Jdata(N_multi+1+n)=-Jdata(N_multi+1+n+2)+(2*n+3)/z*Jdata(N_multi+1+n+1);
% end
% 
% end
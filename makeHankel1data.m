function Hdata = makeHankel1data(N_multi,z)
% Gives the spherical Hankel function h_l^(1) evaluated at z. The first two
% (l=0, l=1) are defined explicitly then a recurrence relation is used to
% generate the remaining l=-N_multi,...,N_multi

Hdata = zeros(N_multi+1,1);

Hdata(1) = sqrt(pi/2/z)*besselh(1/2,1,z);
if N_multi > 0
    Hdata(2) = sqrt(pi/2/z)*besselh(3/2,1,z);

    for n=3:N_multi+1
        Hdata(n)=-Hdata(n-2)+(2*n-3)/z*Hdata(n-1);
    end
end

end


% Hdata = zeros(2*N_multi+1,1);
% 
% Hdata(N_multi+1+0) = sqrt(pi/2/z)*besselh(1/2,1,z);
% Hdata(N_multi+1+1) = sqrt(pi/2/z)*besselh(3/2,1,z);
% 
% for n=2:N_multi
%     Hdata(N_multi+1+n)=-Hdata(N_multi+1+n-2)+(2*n-1)/z*Hdata(N_multi+1+n-1);
% end
% 
% for n=-1:-1:-N_multi
%     Hdata(N_multi+1+n)=-Hdata(N_multi+1+n+2)+(2*n+3)/z*Hdata(N_multi+1+n+1);
% end
% 
% end
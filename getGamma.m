function [Res, DS1Res, DS2Res, DR1Res, DR2Res, ...
    DS1DR1Res, DS1DR2Res, DS2DR1Res, DS2DR2Res, IsInfIdx] = getGamma(Source, Receiver, k)
%Pre:  Source is a (2 x N) array, Source is not part of the boundary,
%      Receiver is a (2 x M) array, Receiver is not part of the boundary,
%      k is a real positive number
%Post: res is (M x N) array. Same for the others
%      IsInfIdx is a column of positive integers
%
% res represents the Gamma Function at the N source points.
% IsInfIdx show almost inf entries in res, i.e. res(IsInfIdx)) is large
% Source is not considered to be inside the domain.

twopi=2*pi;
eg=0.57721566490153286060651209008240243104215933593992;

M=size(Receiver,2);
N=size(Source,2);

SmR=(Source(1,:).'-Receiver(1,:)).^2 +(Source(2,:).'-Receiver(2,:)).^2;
IsInfIdx=find(SmR < 1e-14);

Res = -1i/4*besselh(0,k*sqrt(SmR.'));
logSmRplus=log(sqrt(SmR(mod(IsInfIdx,M*N)+1)));
logSmRminus=log(sqrt(SmR(mod(IsInfIdx-2,M*N)+1)));
logSmRplus(isinf(logSmRplus))=logSmRminus(isinf(logSmRplus));
logSmRminus(isinf(logSmRminus))=logSmRplus(isinf(logSmRminus));
logSmR = (logSmRplus+logSmRminus)/2;
Res(IsInfIdx) = (logSmR-2)/twopi+((log(k/2)+eg)/twopi-1i/4);
% We assign to inf number values which allows to integrate over a
% logarithm integral. Only thought about the case in which Source ==
% Receiver

%From here on we assume that sources and receivers are not equal
SmR=SmR.';
DS1Res = 1i*k/4*besselh(1,k*sqrt(SmR)).*(Source(1,:) - Receiver(1,:).')./sqrt(SmR);
DS2Res = 1i*k/4*besselh(1,k*sqrt(SmR)).*(Source(2,:) - Receiver(2,:).')./sqrt(SmR);
DR1Res = 1i*k/4*besselh(1,k*sqrt(SmR)).*(Receiver(1,:) .'- Source(1,:))./sqrt(SmR);
DR2Res = 1i*k/4*besselh(1,k*sqrt(SmR)).*(Receiver(2,:).' - Source(2,:))./sqrt(SmR);

DS1DR1Res = -1i/4*(k*(Source(1,:)-Receiver(1,:).').^2.*besselh(0,k*sqrt(SmR))./SmR...
    -((Source(1,:)-Receiver(1,:).')+(Source(2,:)-Receiver(2,:).')).*((Source(1,:)-Receiver(1,:).')-(Source(2,:)-Receiver(2,:).')).*besselh(1,k*sqrt(SmR))./sqrt(SmR).^3);
DS1DR2Res = 1i*k/4.*(Source(1,:) - Receiver(1,:).').*(Source(2,:) - Receiver(2,:).').*besselh(2,k*sqrt(SmR))./SmR;
DS2DR1Res = 1i*k/4.*(Source(1,:) - Receiver(1,:).').*(Source(2,:) - Receiver(2,:).').*besselh(2,k*sqrt(SmR))./SmR;
DS2DR2Res = -1i/4*(k*(Source(2,:)-Receiver(2,:).').^2.*besselh(0,k*sqrt(SmR))./SmR...
    +((Source(1,:)-Receiver(1,:).')+(Source(2,:)-Receiver(2,:).')).*((Source(1,:)-Receiver(1,:).')-(Source(2,:)-Receiver(2,:).')).*besselh(1,k*sqrt(SmR))./sqrt(SmR).^3);
end

%TEST VECTORISATION
% A=zeros(2,3);
% Dom=shape.Ellipse(0.1,0.1,1000)+[0;0.3];
% Source=[0.1,-0.1;0.1,0.1];
% Receiver=[0.1,-0.1, 0.05;0.5,0.5, 0.6];
% 
% AFull=getNeumannFct(Dom, Source, Receiver, k);
% k=1.5;
% for i=1:size(A,1)
%     for j=1:size(A,2)
%         A(i,j) = getNeumannFct(Dom, Source(:,i), Receiver(:,j), k);
%     end
% end
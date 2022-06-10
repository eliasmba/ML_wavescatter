function [y, nor]= optimalDiskLoc(U, Y, z, x, k, N)
%PRE: U is a n x size(z,2) real matrix
%     Y is a n x size(x,2) real positive matrix
%     z is a 2 x mz coordinates with points in (-inf,-1.0) x R
%     x is a 2 x mx coordinates with points in (1.0, inf) x R
%     k is a real positive number
%     N is a positive integer
%Post: y are N points in (-1,1) x (-1,1)
%Desc: Uses Foldy-Lax approximation to determine points in (-1,1) x (-1,1)
%such that  ||U*N-Y|| is minimized.

YtUT = (pinv(Y)*U).';

absDy = @(y)sqrt((y(1,:).'-y(1,:)).^2+(y(2,:).'-y(2,:)).^2);
absDx = @(y)sqrt((y(1,:).'-x(1,:)).^2+(y(2,:).'-x(2,:)).^2);
absDz = @(y)sqrt((y(1,:).'-z(1,:)).^2+(y(2,:).'-z(2,:)).^2);
Gamma = @(y)-1i/4*besselh(0,k*absDy(y));
Gammax = @(y)-1i/4*besselh(0,k*absDx(y));
Gammaz = @(y)-1i/4*besselh(0,k*absDz(y));
D1Gamma = @(y)1i*k/4*besselh(1,k*absDy(y)).*(y(1,:).'-y(1,:))./absDy(y);
D2Gamma = @(y)1i*k/4*besselh(1,k*absDy(y)).*(y(2,:).'-y(2,:))./absDy(y);
D1Gammax = @(y)1i*k/4*besselh(1,k*absDx(y)).*(y(1,:).'-x(1,:))./absDx(y);
D1Gammaz = @(y)1i*k/4*besselh(1,k*absDz(y)).*(y(1,:).'-z(1,:))./absDz(y);
D2Gammax = @(y)1i*k/4*besselh(1,k*absDx(y)).*(y(2,:).'-x(2,:))./absDx(y);
D2Gammaz = @(y)1i*k/4*besselh(1,k*absDz(y)).*(y(2,:).'-z(2,:))./absDz(y);

%f = @(y)N*eye(N)+Gamma(y)+Gammax(y).'*YtU*Gammaz(y);
%F = @(y)sum(f(y).^2,'all'); %y= [0.2, 0.3, -0.5, ...;
                             %  -0.1, 0.4, -0.1, ...]

% DF = @(y)[2*( (D1Gamma(y)(i,:)+D1Gammax(y)(i,:).'*YtU*Gammaz(y)).*f(y)(i,:).'...
%             + f(y)(:,i).'.*(D1Gamma(y)(:,i)+Gammax(y).'*YtU*D1Gammaz(y)(:,i))  );...
%           2*( (D1Gamma(y)(i,:)+D1Gammax(y)(i,:).'*YtU*Gammaz(y)).*f(y)(i,:).'...
%             + f(y)(:,i).'.*(D1Gamma(y)(:,i)+Gammax(y).'*YtU*D1Gammaz(y)(:,i))  ) ];

% Initial Guess
y0 = (rand(2,N)-0.5)*2;

% Gradient descent
y = y0;
gamma = 3;
for i = 1:20 
    GammaMat = Gamma(y);
    GammaxMat = Gammax(y);
    GammazMat = Gammaz(y);
    D1GammaMat = D1Gamma(y);
    D2GammaMat = D2Gamma(y);
    D1GammaxMat = D1Gammax(y);
    D1GammazMat = D1Gammaz(y);
    D2GammaxMat = D2Gammax(y);
    D2GammazMat = D2Gammaz(y);
    
    GammaMat(1:1+size(GammaMat,1):end) = 0;
%     GammaxMat(1:1+size(GammaxMat,1):end) = 0;
%     GammazMat(1:1+size(GammazMat,1):end) = 0;
    D1GammaMat(1:1+size(D1GammaMat,1):end) = 0;
    D2GammaMat(1:1+size(D2GammaMat,1):end) = 0;
%     D1GammaxMat(1:1+size(D1GammaxMat,1):end) = 0;
%     D1GammazMat(1:1+size(D1GammazMat,1):end) = 0;
%     D2GammaxMat(1:1+size(D2GammaxMat,1):end) = 0;
%     D2GammazMat(1:1+size(D2GammazMat,1):end) = 0;
    
    RefMat = real(eye(N)+GammaMat/N+GammazMat*YtUT*GammaxMat.'/N);
    ImfMat = imag(eye(N)+GammaMat/N+GammazMat*YtUT*GammaxMat.'/N);
    
    DFMat = 2*[(reshape(sum(RefMat.*real(D1GammaMat+D1GammazMat*YtUT*GammaxMat.'),2),1,N)/N...
        +reshape(sum(RefMat.*real(D1GammaMat+GammazMat*YtUT*D1GammaxMat.'),1),1,N)/N...
        +reshape(sum(ImfMat.*imag(D1GammaMat+D1GammazMat*YtUT*GammaxMat.'),2),1,N)/N...
        +reshape(sum(ImfMat.*imag(D1GammaMat+GammazMat*YtUT*D1GammaxMat.'),1),1,N)/N);
        (reshape(sum(RefMat.*real(D2GammaMat+D2GammazMat*YtUT*GammaxMat.'),2),1,N)/N...
        +reshape(sum(RefMat.*real(D2GammaMat+GammazMat*YtUT*D2GammaxMat.'),1),1,N)/N...
        +reshape(sum(ImfMat.*imag(D2GammaMat+D2GammazMat*YtUT*GammaxMat.'),2),1,N)/N...
        +reshape(sum(ImfMat.*imag(D2GammaMat+GammazMat*YtUT*D2GammaxMat.'),1),1,N)/N)];
%     DFMat = 2*[(reshape(sum(RefMat.*real(D1GammaMat+D1GammazMat*YtUT*GammaxMat.'),2),1,N)/N...
%         +reshape(sum(RefMat.*real(D1GammaMat+GammazMat*YtUT*D1GammaxMat.'),1),1,N)/N);
%         (reshape(sum(RefMat.*real(D2GammaMat+D2GammazMat*YtUT*GammaxMat.'),2),1,N)/N...
%         +reshape(sum(RefMat.*real(D2GammaMat+GammazMat*YtUT*D2GammaxMat.'),1),1,N)/N)];
           
    %disp(num2str(sum((RefMat).^2+(ImfMat).^2,'all')));
    %disp(num2str(sum((RefMat).^2,'all')));
    y = y - gamma*DFMat;
end
nor=sum((RefMat).^2+(ImfMat).^2,'all');%
%nor=sum((RefMat).^2,'all');
%disp(num2str(nor));
end


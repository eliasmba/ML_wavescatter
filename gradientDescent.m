function[y, nor] = gradientDescent(U, Y, z, x,R,kappa0, kappa_b, rho_b, rho0, k, N)

YtUT = (pinv(Y)*U).';

% absDy = @(y)sqrt((y(1,:).'-y(1,:)).^2+(y(2,:).'-y(2,:)).^2);
% absDx = @(y)sqrt((y(1,:).'-x(1,:)).^2+(y(2,:).'-x(2,:)).^2);
% absDz = @(y)sqrt((y(1,:).'-z(1,:)).^2+(y(2,:).'-z(2,:)).^2
% Gamma = @(y)-1i/4*besselh(0,k*absDy(y));
% Gammax = @(y)-1i/4*besselh(0,k*absDx(y));
% Gammaz = @(y)-1i/4*besselh(0,k*absDz(y));
% D1Gamma = @(y)1i*k/4*besselh(1,k*absDy(y)).*(y(1,:).'-y(1,:))./absDy(y);
% D2Gamma = @(y)1i*k/4*besselh(1,k*absDy(y)).*(y(2,:).'-y(2,:))./absDy(y);
% D1Gammax = @(y)1i*k/4*besselh(1,k*absDx(y)).*(y(1,:).'-x(1,:))./absDx(y);
% D1Gammaz = @(y)1i*k/4*besselh(1,k*absDz(y)).*(y(1,:).'-z(1,:))./absDz(y);
% D2Gammax = @(y)1i*k/4*besselh(1,k*absDx(y)).*(y(2,:).'-x(2,:))./absDx(y);
% D2Gammaz = @(y)1i*k/4*besselh(1,k*absDz(y)).*(y(2,:).'-z(2,:))./absDz(y);

y0 = (rand(3,N)-0.5)*2;

% Gradient descent
y = y0;
gamma = 0.0000000000000000000001;
for i = 1:3 
    disp(i);
    disp('Here');
    disp(y);
%     GammaMat = Gamma(y);
%     GammaxMat = Gammax(y);
%     GammazMat = Gammaz(y);
%     D1GammaMat = D1Gamma(y);
%     D2GammaMat = D2Gamma(y);
%     D1GammaxMat = D1Gammax(y);
%     D1GammazMat = D1Gammaz(y);
%     D2GammaxMat = D2Gammax(y);
%     D2GammazMat = D2Gammaz(y);

 
    
%     GammaMat(1:1+size(GammaMat,1):end) = 0;
%     D1GammaMat(1:1+size(D1GammaMat,1):end) = 0;
%     D2GammaMat(1:1+size(D2GammaMat,1):end) = 0;

    [D1GammaMat, D2GammaMat, D3GammaMat, D1GammaxMat, D2GammaxMat, D3GammaxMat, D1GammazMat, D2GammazMat, D3GammazMat] = derivatives(z, x, y, U, R*ones(size(y,2),1), 0, kappa0, rho0, kappa_b, rho_b);
    
    [GammaMat, GammaxMat, GammazMat] = get_functions(z, x, y, U, R*ones(size(y,2),1), 0, kappa0, rho0, kappa_b, rho_b);
    GammaMat(1:1+size(GammaMat,1):end) = 0;
    
%     disp(size(YtUT));
%     disp(size(GammazMat));
    
    RefMat = real(eye(N)+GammaMat.'/N+GammazMat*YtUT*GammaxMat.'/N);
    ImfMat = imag(eye(N)+GammaMat.'/N+GammazMat*YtUT*GammaxMat.'/N);
    
    DFMat = 2*[(reshape(sum(RefMat.*real(D1GammaMat+D1GammazMat*YtUT*GammaxMat.'),2),1,N)/N...
        +reshape(sum(RefMat.*real(D1GammaMat+GammazMat*YtUT*D1GammaxMat.'),1),1,N)/N...
        +reshape(sum(ImfMat.*imag(D1GammaMat+D1GammazMat*YtUT*GammaxMat.'),2),1,N)/N...
        +reshape(sum(ImfMat.*imag(D1GammaMat+GammazMat*YtUT*D1GammaxMat.'),1),1,N)/N);
        (reshape(sum(RefMat.*real(D2GammaMat+D2GammazMat*YtUT*GammaxMat.'),2),1,N)/N...
        +reshape(sum(RefMat.*real(D2GammaMat+GammazMat*YtUT*D2GammaxMat.'),1),1,N)/N...
        +reshape(sum(ImfMat.*imag(D2GammaMat+D2GammazMat*YtUT*GammaxMat.'),2),1,N)/N...
        +reshape(sum(ImfMat.*imag(D2GammaMat+GammazMat*YtUT*D2GammaxMat.'),1),1,N)/N);
        (reshape(sum(RefMat.*real(D3GammaMat+D3GammazMat*YtUT*GammaxMat.'),2),1,N)/N...
        +reshape(sum(RefMat.*real(D3GammaMat+GammazMat*YtUT*D3GammaxMat.'),1),1,N)/N...
        +reshape(sum(ImfMat.*imag(D3GammaMat+D3GammazMat*YtUT*GammaxMat.'),2),1,N)/N...
        +reshape(sum(ImfMat.*imag(D3GammaMat+GammazMat*YtUT*D3GammaxMat.'),1),1,N)/N)];
    
    y = y - gamma*DFMat;
end
nor=sum((RefMat).^2+(ImfMat).^2,'all');
end
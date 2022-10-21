clear
format compact

%% MNIIST WORKING SET
[XTrain,YTrain,anglesTrain] = digitTrain4DArrayData;
[XTest,YTest,anglesTest] = digitTest4DArrayData;

disp(size(XTrain))
disp(size(YTrain))

disp(size(XTest))

Y = zeros(5000,10);
Ytest = zeros(5000,10);
labels = 0:9;
for n=1:5000
    Y(n,labels(YTrain(n))+1) = 1;
    Ytest(n,labels(YTest(n))+1) = 1;
end

U = reshape(XTrain,28*28,5000).';
Util = reshape(XTest,28*28,5000).';

%% Reducing
labels = 0:2;
colYidx = Y(:,1)==1 | Y(:,2)==1 | Y(:,3)==1;
U = U(colYidx,:);
Y = Y(colYidx,1:3);
colYtestidx = Ytest(:,1)==1 | Ytest(:,2)==1 | Ytest(:,3)==1;
Util = Util(colYtestidx,:);
Ytest = Ytest(colYtestidx,1:3);

%randomize test set
randidx = randperm(size(Ytest,1));
Util = Util(randidx,:);
Ytest = Ytest(randidx,:);
N = 10;
% Define size of resonators
R = 0.1;
z = [-1.1*ones(1,size(U,2)); linspace(-1.0,1.0,size(U,2)); 0.5*ones(1, size(U, 2))];
x = [+1.1*ones(1,size(Y,2)); linspace(-0.7,0.7,size(Y,2));0.5*ones(1, size(Y,2))];

%%% Material parameters
high = 5000;

rho0 = high;             % density of background
kappa0 = high;           % bulk modulus of background
v = sqrt(kappa0/rho0);  % speed of sound in background

rho_b = 1;            % density of resonators  
kappa_b = 1;          % bulk modulus of resonators
v_b = sqrt(kappa_b/rho_b);  % speed of sound in air

% High contrast parameters \delta
delta=rho_b/rho0;
omega=3;
y = (rand(3,3)-0.5)*2; %initial guess
[M, Sigmax, Sigmaz] = get_functions(z, x, y, U, R*ones(size(y,2),1), 0, kappa0, rho0, kappa_b, rho_b);


%% Test Outcome 
k0=omega*sqrt(rho0/kappa0);
kb=omega*sqrt(rho_b/kappa_b);
N=10;
disp(['k0 : ', num2str(k0), 'kb : ', num2str(kb)]);
%y=rand(2,N);
[y, nor] = gradientDescent(U, Y, z, x, rho_b, rho0, omega, N);


Nk = Sigmax(y).'*(M\Sigmaz(y));
norUY = sum(abs(U*Nk.'-Y), 'all'); 
disp(['nor = ', num2str(norUY)]);
plot(y(1,:),y(2,:),'.g','MarkerSize',5); hold on;
plot(z(1,:),z(2,:),'.r','MarkerSize',5); 
plot(x(1,:),x(2,:),'.b','MarkerSize',5); 
axis([-1.2, 1.2, -1.2, 1.2]);hold off;

% Applying Test Set

accuracy = zeros(1,length(labels)); %preload

Ytil = abs(Util*Nk.');
sigCostTemp = 1./(1+exp(-abs(Ytil)));
[~,estIdxLabel] = max((Ytil).');
estLabel = (estIdxLabel-1).';
for q=labels
    qIdx = (Ytest(:,q+1) == 1);
    accuracy(q+1) = sum(estLabel(qIdx) == q)/sum(qIdx);
    disp(['Accuracy of label ',num2str(q), ' is ', num2str(accuracy(q+1))]);
end
disp(['Mean Accuracy : ', num2str(mean(accuracy))]);
disp('end');

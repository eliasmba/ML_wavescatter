clear, close all

L = 35e-3;              % length of cochlea
N = 10;                 % number of resonators
L_end = L - 0.005;

s = 1.05;
Rad = L_end*(1-s)/(1-s^N)/3;
for i = 1:N
    R(i) = Rad*s^(i-1);
end

%%% Material parameters
rho0 = 1e3;             % density of water
kappa0 = 2e9;           % bulk modulus of water
v = sqrt(kappa0/rho0);  % speed of sound in water

rho_b = 1.2;            % density of resonators  
kappa_b = 1e5;          % bulk modulus of resonators
v_b = sqrt(kappa_b/rho_b);  % speed of sound in air


% High contrast parameters \delta
delta=rho_b/rho0;

% Define positions of resonators
cx = R(1)*ones(1,N);
if N > 1
    for i = 2:N
        cx(i) = cx(i-1) + 2*R(i-1) + R(i);
    end
end

%%% Plot the geometry
cy = zeros(1,N); cz = zeros(1,N);
figure, hold on
t = linspace(0,2*pi);
for n = 1:length(R)
    plot(cx(n)+R(n)*cos(t), cy(n)+R(n)*sin(t),'k')
end
daspect([1 1 1])
hold off

%% Compute the resonances using the capacitance matrix, which is approximated using the multipole expansion method
C = capacitance(cx,cy,cz,R,rho0,rho_b,kappa0,kappa_b,delta);
size(C)
Cvol = 3/4/pi*diag(R.^-3)*C;
capres = sqrt(delta*v_b^2*eig(real(Cvol)));

%% Plot the resonances
figure
hold on
scatter(capres,ones(1,N),'ob')
ylim([0.85, 1.2]); set(gca,'ytick',[]);
legend('full multipole method')

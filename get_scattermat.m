function scattermat = get_scattermat(m, p, n, k, k0, N_Multi)
R = 0.015;
%% Build Up Geometry
% Note: While Loops might not end
%Disks
nd_loop = 1;
centers = zeros(3,n);
while nd_loop <= n
    newCoord = (rand(3,1)-0.5).*2; % Vector in [-1.0,1.0]^2
    if ~any(sum((centers(:,1:nd_loop-1)-newCoord).^2)<(R*sqrt(2)*2+Disk2DiskThresh)^2) % if newCoord are close to an existing disk or the origin
        centers(:,nd_loop) = newCoord;
        nd_loop=nd_loop+1;
    end
end

% Sources
nr_loop = 1;
source = zeros(3,p);
while nr_loop<=p
    zrnd2 = (rand(2,1)-0.5)*2;
    zrnd1 = rand(1,1)*0.9-1;
    if ~any(sum((centers-[zrnd1;zrnd2]).^2)<(R*sqrt(2)+Disk2DiskThresh/2)^2) ...
      && ~any(sum((sources(:,1:nr_loop-1)-[zrnd1;zrnd2]).^2)<Pt2PtThresh^2)
        sources(:,nr_loop) = [zrnd1;zrnd2];
        nr_loop = nr_loop + 1;
    end
end

%Receivers
nr_loop = 1;
receivers= zeros(3,m);
while nr_loop<=m
    xrnd2 = (rand(2,1)-0.5)*2;
    xrnd1 = (rand(1,1)*0.9+0.1);
    if  ~any(sum((centers-[xrnd1;xrnd2]).^2)<(R*sqrt(2)+Disk2DiskThresh/2)^2) ...
      && ~any(sum((receivers(:,1:nr_loop-1)-[xrnd1;xrnd2]).^2)<Pt2PtThresh^2)
        receivers(:,nr_loop) = [xrnd1;xrnd2];
        nr_loop = nr_loop + 1;
    end
end

%% evaluate scattering
scattermat = zeros(p, m);
for i=1:p
    scattermat(i,:)= LFRsolver(sources, receivers, eye(i), centers, R, k, k0, N_Multi);
end
end

% Clean memory and command window
clear,clc,close all

%% Parameters
N = 100;             % Number of nodes
W = 200;             % length of the network
L = 200;             % width of the network
Ei = 2;              % Initial energy of each node (joules)
CHpl = 3000;         % Packet size for cluster head per round (bits)

p = 5/100;           % desired percentage of cluster heads
R = 50;              % Range for cluster
pMin = 10^-4;        % Lowest possible CH_prop

num_rounds = 2000;   % Max Number of simulated rounds
NonCHpl = 200;       % Packet size for normal node per round (bits)
Tsetup = 4;          % average Time in seconds taken in setup phase
Tss = 10;            % average Time in seconds taken in steady state phase
Etrans = 1.0000e-05; % Energy for transmitting one bit 
Erec = 1.0000e-05;   % Energy for receiving one bit 
Eagg = 1.0000e-07;   % Data aggregation energy
Efs = 0.3400e-9;     % Energy of free space model amplifier
% Position of sink
SX = W/2; SY = L/2;


%% First Leach algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%

% 1st row: states of being a CH, 1:never been CH, 0:has been CH 
% 2nd: x-position, 3rd: y-position
net = [ones(1,N);rand([1,N])*W;rand([1,N])*L];

% Preallocation for energy calculations
E = Ei*ones(1,N);          % Energy left in each node
EH = zeros(1,num_rounds); 
% Preallocation for dead nodes calculations
Ecrit = 0;                 % Critical energy left in node to call it alive
Ncrit = fix((95/100)*N);   % Critical number for dead nodes to stop simulation
Dead = false(1,N);         % Status of all nodes 0:Alive 1:Dead
DeadH = zeros(1,num_rounds);
% Preallocation for Bits sent calculations
BitsH = zeros(1,num_rounds);
figure('Position',[34 30 792 613]);

% Simulating for each round
for r=1:num_rounds % iterating on each round
    
    %%%% Choosing Clusters heads %%%%
    [net(1,:),CH] = Leach_algo(net(1,:),Dead,p,r);
    tmp = find(CH);
    for i=1:N
        if isempty(net(2,CH))
            
        else
            [~,aa]=min(sqrt((net(2,CH) - net(2,i)).^2 + (net(3,CH) - net(3,i)).^2));
            net(1,i) = tmp(aa);
        end
    end
    %%%% Energy calculations %%%%
    EH(r) = sum(E); %get all energy left in all nodes
    % first CH
    numClust = length(find(CH));
    D = sqrt((net(2,CH) - SX).^2 + (net(3,CH) - SY).^2);
    E(CH) = E(CH) - (((Etrans+Eagg)*CHpl)+(Efs*CHpl*(D.^ 2))+(NonCHpl*Erec*round(N/numClust)));
    % second rest of nodes
    rest = N-numClust-sum(double(Dead));
    mD = zeros(1,rest); tmp = net(2:3,~CH&~Dead);
    for i=1:rest, mD(i) = fun(tmp(1,i),tmp(2,i),net,CH,SX,SY); end
    E(~CH&~Dead) = E(~CH&~Dead) - ((NonCHpl*Etrans) + (Efs*CHpl*(mD.^2)) + ((Erec+Eagg)*CHpl));
    %finally updating alive status to all nodes
    E(Dead) = 0;
    Dead(E<=Ecrit) = true ; DeadH(r)=sum(double(Dead));
    %%%% sent bits %%%%
    BitsH(r+1) = BitsH(r) + numClust*CHpl + rest*NonCHpl;
    
    %%%% Showing updated net %%%%
    net = DrawNet(net,N,CH,Dead,SX,SY,1);
    title(['Normal nodes:Black ---- CH:Red ---- Dead:Empty circle --- round (',num2str(r),')']);
    drawnow
    
    if DeadH(r)>=Ncrit,break;end % Stop simulation when 5% or less is alive
end
close all
T_L = (Tsetup+Tss)*(0:r-1);
EH = EH(1:r); EHdis_L = (N*Ei)-EH;
DeadH = DeadH(1:r); AliveH_L = N-DeadH;
BitsH_L = BitsH(2:r+1);

%% Second HEED algorithm %%%%%%%%%%%%%%%%%%%%%%%

% 1st row: Clustering indexing 
% 2nd: x-position, 3rd: y-position
net = [zeros(1,N);net(2:3,:)];

% calculating costs
cost = zeros(1,N);
for i=1:N
    Dist = sqrt(((net(2,:)-net(2,i)).^2) + ((net(3,:)-net(3,i)).^2));
    Snbr = Dist <= R;
    cost(i) = sum(Dist(Snbr))/(sum(Snbr)-1);
end

% Preallocation for energy calculations
E = Ei*ones(1,N);          % Energy left in each node
EH = zeros(1,num_rounds); 
% Preallocation for dead nodes calculations
Ecrit = 0;                 % Critical energy left in node to call it alive
Ncrit = fix((98/100)*N);   % Critical number for dead nodes to stop simulation
Dead = false(1,N);         % Status of all nodes 0:Alive 1:Dead
DeadH = zeros(1,num_rounds);
% Preallocation for Bits sent calculations
BitsH = zeros(1,num_rounds);
figure('Position',[34 30 792 613]);

% Simulating for each round
for r=1:num_rounds % iterating on each round
    
    %%%% Choosing Clusters heads %%%%
    [CH,net] = HEED_algo(R,Dead,p,pMin,E,Ei,net,cost);
    
    %%%% Energy calculations %%%%
    EH(r) = sum(E); %get all energy left in all nodes
    % first CH
    numClust = length(find(CH));
    D = sqrt((net(2,CH) - SX).^2 + (net(3,CH) - SY).^2);
    E(CH) = E(CH) - (((Etrans+Eagg)*CHpl)+(Efs*CHpl*(D.^ 2))+(NonCHpl*Erec*round(N/numClust)));
    % second rest of nodes
    rest = N-numClust-sum(double(Dead));
    mD = zeros(1,rest); tmp = net(2:3,~CH&~Dead);
    for i=1:rest, mD(i) = funH(tmp(1,i),tmp(2,i),net,CH,SX,SY); end
    E(~CH&~Dead) = E(~CH&~Dead) - ((NonCHpl*Etrans) + (Efs*CHpl*(mD.^2)) + ((Erec+Eagg)*CHpl));
    %finally updating alive status to all nodes
    E(Dead) = 0; CH(Dead) = false;
    Dead(E<=Ecrit) = true ; DeadH(r)=sum(double(Dead));
    %%%% sent bits %%%%
    BitsH(r+1) = BitsH(r) + numClust*CHpl + rest*NonCHpl;
    
    %%%% Showing updated net %%%%
    net = DrawNet(net,N,CH,Dead,SX,SY,1);
    title(['Normal nodes:Black ---- CH:Red ---- Dead:Empty circle --- round (',num2str(r),')']);
    drawnow
    if sum(Dead)>=Ncrit,break;end % Stop simulation when 5% or less is alive
end
T = (Tsetup+Tss)*(0:r-1);
EH = EH(1:r); EHdis = (N*Ei)-EH;
DeadH = DeadH(1:r); AliveH = N-DeadH;
BitsH = BitsH(2:r+1);
close all

%% Plotting analysis of network preformance & comparison

figure('Position',[131 59 792 613]);
plot(T,EHdis,'-x')
plot(T,EHdis,'-x',T_L,EHdis_L,'-x');
xlabel('Time (s)'); ylabel('Energy (j)')
title('Total energy dissipated')
legend('HEED','LEACH')

figure('Position',[298 66 792 613]);
plot(T,AliveH,'-x')
plot(T,AliveH,'-x',T_L,AliveH_L,'-x'); 
xlabel('Time (s)'); ylabel('No of live Sensors nodes')
title('Life time of sensor nodes')
legend('HEED','LEACH')

figure('Position',[511 61 792 613]);
plot(T,BitsH,'-x')
plot(T,BitsH,'-x',T_L,BitsH_L,'-x'); 
xlabel('Time (s)'); ylabel('Throughput (no of packets)')
title(['Throughput (' num2str(N) ' Sensor nodes)'])
legend('HEED','LEACH')


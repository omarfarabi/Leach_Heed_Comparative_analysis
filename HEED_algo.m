function [CH_F,net] = HEED_algo(R,D,p,pMin,E,Emax,net,cost)

% R : range of a cluster
% D : State of nodes (Dead or alive)
% p : Cprop
% pMin : Minimum value of CHprop
% E    : Residual energy in all nodes
% Emax : Maximum energy each node can hold
% net  : Matrix holding position of each node and their clustering
% cost : cost for each node to ba a cluster head

% CH : vector showing the final Cluster Heads
N = size(net,2); % number of nodes

% initializing
CH_prop = max([p*(E./Emax);pMin*ones(1,N)]);
clust_idx = zeros(1,N); % clusters indices
CH_F = false(1,N);    % Final CH
CH_tent = false(1,N); % Tentative CH

% Clustering
for i=1:N
    if ~D(i) % to make sure node is not dead
        % Clustering rang of this node
        Dist = sqrt(((net(2,:)-net(2,i)).^2) + ((net(3,:)-net(3,i)).^2));
        Snbr = Dist <= R;
        CH_prev = 0;
        
        % Repeat
        while CH_prev~=1
            tmp = rand;
            if sum((CH_tent&Snbr&(~D)))>0
                tmp = cost; tmp(~(CH_tent&Snbr&(~D)))= inf;
                [~,my_CH] = min(tmp);
                clust_idx(i) = my_CH;
                if my_CH==i
                    if CH_prop(i)==1
                        CH_F(i)=true;
                    else
                        CH_tent(i)=true;
                    end
                end
            elseif CH_prop(i)==1
                CH_F(i)=true;
                clust_idx(i)=i;
            elseif tmp <= CH_prop(i)
                CH_tent(i)=true;
            end
            CH_prev = CH_prop(i);
            CH_prop(i) = min(2*CH_prop(i),1);
        end
        
        % finalize
        if ~CH_F(i)
            if sum((CH_tent&Snbr&(~D)))>0
                tmp = cost; tmp(~(CH_tent&Snbr&(~D)))= inf;
                [~,my_CH] = min(tmp);
                clust_idx(i) = my_CH;
            else
                CH_tent(i)=true;
                CH_F(i)=true;
                clust_idx(i)=i;
            end
        end
    end
end

net(1,:) = clust_idx;
CH_F(D) = false;







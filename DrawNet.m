function net = DrawNet(net,N,CH,D,SX,SY,ShowClust)

% This function plots the shape of the network
% showing all nodes and its status
% If there is cluster map [(3-D matrix)
% where first and second dimensions are
% x and y postions of cluster border
% and the third dimension is number of clusters]
% then the function plots those clusters in green lines

% normal nodes
plot(net(2,~CH&~D),net(3,~CH&~D),'ko','MarkerSize',...
    5,'MarkerFaceColor','k'); hold on;
% dead nodes
plot(net(2,D),net(3,D),'ko','MarkerSize',5);
% Cluster heads
plot(net(2,CH),net(3,CH),'ko','MarkerSize',...
    7,'MarkerFaceColor','r');
% The sink
scatter(SX,SY,100,'bo','filled')
text(SX+2,SY+2,'Sink','FontSize',10,'VerticalAlignment','Baseline');
% titles
s = int2str((1:N)');
text(net(2,:)+1,net(3,:)+1,s,'FontSize',8,'VerticalAlignment','Baseline');
xlabel('\it x \rm [m] \rightarrow');
ylabel('\it y \rm [m] \rightarrow');

if ShowClust
    idx = net(1,:); Clust = unique(idx); Clust(Clust==0)=[];
    X = net(2:3,:); X = X';
    for i=1:length(Clust)
        tmp= idx==Clust(i); tmp=X(tmp,:);
        for j=1:size(tmp,1)
            plot([tmp(j,1),X(Clust(i),1)],[tmp(j,2),X(Clust(i),2)],'-','Color','g')
        end
    end
end

hold off
function d = funH(x,y,net,CH,SX,SY)
if isempty(net(2,CH))
    d=sqrt((x - SX).^2 + (y - SY).^2);
else
    idx = (net(2,:)==x) & (net(3,:)==y);
    cluster = net(1,:) == net(1,idx);
    idx = cluster & CH;
    
    d=sqrt((x - net(2,idx)).^2 + (y - net(3,idx)).^2);
end
function d = fun(x,y,net,CH,SX,SY)
if isempty(net(2,CH))
    d=sqrt((x - SX).^2 + (y - SY).^2);
else
    d=min(sqrt((net(2,CH) - x).^2 + (net(3,CH) - y).^2));
end
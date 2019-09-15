function [G,CH] = Leach_algo(G,D,p,r)

tmp = rand(size(G)); %random number for each node
T = G; 
idx = T==1;
if isempty(find(idx, 1))
    G = ones(size(G));
    idx = G==1;
end

T(idx) = p / (1-p * mod(r, round(1 / p)));
T(~idx) = 0;
T(D) = 0;
CH = tmp<T; 
G(CH) = 0; 
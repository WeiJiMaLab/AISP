function [X0,LB,UB,PLB,PUB] = get_bads_bounds()

LB = [0,0,-5,0.0001];
UB = [20,50,5,0.5];
PLB = [2,0,-5,0.0001];
PUB = [10,20,5,0.5];

X0 = PLB + rand(size(LB)) .* (PUB-PLB);
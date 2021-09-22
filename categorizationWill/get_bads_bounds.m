function [X0,LB,UB,PLB,PUB] = get_bads_bounds()

%    [ s1,  s2,  s3,  s4, s5, s6, b0,    b1,   lamb]
LB = [-3, -3, -3, -3, -3, -3, -10, -2, eps];
PLB= [-1,-1,-1,-1,-1,-1,-3,-2,0.01];
PUB= [5,4,3,3,3,2,3,5,0.1];
UB = [6,6,6,6,6,6,10,5,0.25];

X0 = PLB + (PUB - PLB) .* rand(size(PUB));
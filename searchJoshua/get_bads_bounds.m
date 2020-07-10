function [X0,LB,UB,PLB,PUB] = get_bads_bounds()

error('to do')

%    [ s1,  s2,  s3,  s4, b0,    b1,   lamb]
LB = [0.1, 0.1, 0.1, 0.1, -10, 0.01, 0.0001];
PLB= [  1,   1,   1,   1,  -5,  0.1, 0.01];
PUB= [ 25,  25,  25,  25,   5,    5, 0.05];
UB = [100, 100, 100, 100,  10,   25, 0.5];

X0 = PLB + (PUB - PLB) .* rand(size(PUB));
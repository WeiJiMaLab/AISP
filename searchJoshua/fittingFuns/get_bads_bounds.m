function [X0,LB,UB,PLB,PUB] = get_bads_bounds(incSamplesParam)

% INPUT
% incSamplesParam: bool. Whether to include a parameter for the number
%   of samples.

params = {'LnKappa_x', 'LnKappa_x', 'LnKappa_x', 'LnKappa_x', ...
    'Beta0', 'Beta1', 'LapseRate', 'NumSamples'};
LB = [repmat(-6, 1, 4), -10, 0.01, 0.005, 1];
PLB= [repmat(-4, 1, 4),  -5,  0.1, 0.01, 1];
PUB= [repmat(4, 1, 4),   5,    5, 0.2, 1000];
UB = [repmat(7, 1, 4),  10,   25, 1, 1000];

if ~incSamplesParam
    params = params(1: end-1);
    LB = LB(1 : end-1);
    PLB = PLB(1 : end-1);
    PUB = PUB(1 : end-1);
    UB = UB(1 : end-1);
end

X0 = PLB + ((PUB - PLB) .* rand(size(PUB)));

aisp_checkParamOrder(params, incSamplesParam)
    
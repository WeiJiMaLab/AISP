function [X0,LB,UB,PLB,PUB] = get_bads_bounds(incSamplesParam)

% INPUT
% incSamplesParam: bool. Whether to include a parameter for the number
%   of samples.

params = {'LnKappa_x', 'LnKappa_x', 'LnKappa_x', 'LnKappa_x', ...
    'Beta0', 'Beta1', 'LapseRate', 'NumSamples'};
LB = [repmat(-6, 1, 4), -10, 0.01, 0.005, 1];
PLB= [repmat(-4, 1, 4),  -5,  0.1, 0.01, 1];
PUB= [repmat(4, 1, 4),   5,    5, 0.2, 8];
UB = [repmat(7, 1, 4),  10,   25, 1, 12];

if ~incSamplesParam
    params = params(1: end-1);
    LB = LB(1 : end-1);
    PLB = PLB(1 : end-1);
    PUB = PUB(1 : end-1);
    UB = UB(1 : end-1);
end

X0 = PLB + ((PUB - PLB) .* rand(size(PUB)));

% Check that the order the params are specified here matches the ordering used
% for the param vector in the rest of the code
[assumedParamNames, assumedParamOrder] = findParamOrder(incSamplesParam);

for iParam = 1 : length(params)
    
    found = 0;
    for iName = 1 : length(assumedParamNames)
        if ismember(iParam, assumedParamOrder{iName})
            
            found = found +1;
            if ~strcmp(params{iParam}, assumedParamNames{iName})
                error('Bug')
            end
        end
    end
    
    if found ~= 1; error('Bug'); end
end
        
            
    
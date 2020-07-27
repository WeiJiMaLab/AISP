function [X0,LB,UB,PLB,PUB] = get_bads_bounds()

params = {'LnKappa_x', 'LnKappa_x', 'LnKappa_x', 'LnKappa_x', 'Beta0', 'Beta1', 'LapseRate'};
LB = [repmat(-6, 1, 4), -10, 0.01, 0.0001];
PLB= [repmat(-4, 1, 4),  -5,  0.1, 0.01];
PUB= [repmat(4, 1, 4),   5,    5, 0.2];
UB = [repmat(7, 1, 4),  10,   25, 1];

X0 = PLB + ((PUB - PLB) .* rand(size(PUB)));

% Check that the order the params are specified here matches the ordering used
% for the param vector in the rest of the code
[assumedParamNames, assumedParamOrder] = findParamOrder();

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
        
            
    
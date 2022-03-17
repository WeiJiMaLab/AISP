function aisp_checkParamOrder(params, incSamplesParam)
% Check that the order the params as specified here matches the ordering 
% used for the param vector in the rest of the code

% INPUT
% params: cell array of shape 1 x N. N is the number of parameters in 
%   the param vector, and each entry gives the name of the corresponding
%   parameter.
% incSamplesParam: bool. Whether we are including a parameter for the 
%   number of samples.

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
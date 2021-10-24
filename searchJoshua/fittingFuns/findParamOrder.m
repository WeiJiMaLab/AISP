function [paramNames, paramOrder] = findParamOrder(incSamplesParam)
% Sets the order of the parameters in the parameter vector. More precisely,
% paramOrder specifies the columns taken up by a particular parameter in the
% param vector.

% INPUT
% incSamplesParam: bool. Whether to include a parameter for the number
%   of samples.

assert(~isempty(incSamplesParam))

paramNames = {'LnKappa_x', 'Beta0', 'Beta1', 'LapseRate', 'NumSamples'};
paramOrder = {1:4, 5, 6, 7, 8};

if ~incSamplesParam
    paramNames = paramNames(1:end-1);
    paramOrder = paramOrder(1:end-1);
end
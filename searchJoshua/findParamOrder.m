function [paramNames, paramOrder] = findParamOrder()
% Sets the order of the parameters in the parameter vector. More precisely,
% paramOrder specifies the columns taken up by a particular parameter in the
% param vector.

paramNames = {'LnKappa_x', 'Beta0', 'Beta1', 'LapseRate'};
paramOrder = {1:4, 5, 6, 7};

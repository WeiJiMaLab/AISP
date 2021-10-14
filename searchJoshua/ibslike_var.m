function [val_out, var_out, n] = ibslike_var(FUN, PARAMS, RESPMAT, DESIGNMAT, options, var_limit, varargin)
% function ibslike_var(FUN, PARAMS, RESPMAT, DESIGNMAT, options, var_limit)
% This wrapper function runs ibslike with all parameters just passed
% through until the variance of the estimate is smaller than var_limit

% INPUT
% varargin{1}: Boolean. Default false. Whether to display information helpful
% for debuging.

persistent avDuration
persistent numEvals

if length(varargin) >= 1
    debugMode = varargin{1};
else
    debugMode = false;
end

if debugMode
    % Time the duration of the function
    if isempty(avDuration); avDuration = 0; end
    if isempty(numEvals); numEvals = 0; end
    tic
end

[value, var, ~, Details] = ibslike(FUN, PARAMS, RESPMAT, DESIGNMAT, options);

if debugMode
    disp(['Samples per trial used by IBS: ' num2str(Details.NsamplesPerTrial )])
    disp(['Variance from 1 pass of ibslike: ' num2str(var)])
end

val_sum = value;
var_sum = var;
n = 1;

while (var_sum / n / n) > var_limit
    [value, var] = ibslike(FUN, PARAMS, RESPMAT, DESIGNMAT, options);
    var_sum = var_sum + var;
    val_sum = val_sum + value;
    n = n + 1;
end

val_out = val_sum / n;
var_out = var_sum / n / n;

if debugMode
    thisDuration = toc;
    avDuration = ((avDuration*numEvals) + thisDuration) / (numEvals +1);
    numEvals = numEvals +1;
    disp(['Number of ibslike calls: ' num2str(n)])
    disp(['Av time taken by ibslike_var: ' num2str(avDuration) ...
        '. ' num2str(numEvals) ' evals.'])
end


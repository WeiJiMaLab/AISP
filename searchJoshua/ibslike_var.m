function [val_out, var_out, n] = ibslike_var(FUN, PARAMS, RESPMAT, ...
    DESIGNMAT, options, var_limit, varargin)
% function ibslike_var(FUN, PARAMS, RESPMAT, DESIGNMAT, options, var_limit)
% This wrapper function runs ibslike with all parameters just passed
% through until the variance of the estimate is smaller than var_limit

% INPUT
% varargin{1}: Boolean. Default false. Whether to display information helpful
% for debuging.
% varargin{2}: function. A function to call on PARAMS at the very 
% begining of the evaluation. The output will be used as the parameter
% vector instead of PARAMS
% varargin{3}: bool. If true, run a parfor loop, otherwise do not use
%   parallel computation.

persistent avDuration
persistent numEvals

if length(varargin) >= 1
    debugMode = varargin{1};
else
    debugMode = false;
end

if (length(varargin) >= 2) && (~isempty(varargin{2}))
    FunForParams = varargin{2};
    PARAMS = FunForParams(PARAMS);
end

if (length(varargin) >= 3) && (~isempty(varargin{3}))
    runParallel = varargin{3};
else
    runParallel = false;
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
    if ~runParallel
        [value, var] = ibslike(FUN, PARAMS, RESPMAT, DESIGNMAT, options);
        var_sum = var_sum + var;
        val_sum = val_sum + value;
        n = n + 1;
        
    elseif runParallel
        % go to 102% of the samples we expect to need 
        nStep = 1.02 *(var_sum / n / var_limit - n);
        % apply some limits
       	nStep = round(max(min(nStep, 150),1));
        theseVals = nan(nStep, 1);
        theseVars = nan(nStep, 1);
        
        disp(['NUM STEPS: ' num2str(nStep)])
        
        parfor iEval = 1 : nStep
            [thisValue, thisVar] = ibslike(FUN, PARAMS, RESPMAT, ...
                DESIGNMAT, options);
            theseVals(iEval) = thisValue;
            theseVars(iEval) = thisVar;
        end
        
        val_sum = val_sum + sum(theseVals);
        var_sum = var_sum + sum(theseVars);
        n = n + nStep;
    end
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


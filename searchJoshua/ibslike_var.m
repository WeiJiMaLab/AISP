function [val_out, var_out, n] = ibslike_var(FUN, PARAMS, RESPMAT, DESIGNMAT, options, var_limit)
% function ibslike_var(FUN, PARAMS, RESPMAT, DESIGNMAT, options, var_limit)
% This wrapper function runs ibslike with all parameters just passed
% through until the variance of the estimate is smaller than var_limit

[value, var, ~, Details] = ibslike(FUN, PARAMS, RESPMAT, DESIGNMAT, options);
% disp(['Samples per trial: ' num2str(Details.NsamplesPerTrial )])
% disp(sqrt(var))

val_sum = value;
var_sum = var;
n = 1;

while (var_sum / n / n) > var_limit
    [value, var] = ibslike(FUN, PARAMS, RESPMAT, DESIGNMAT, options);
    var_sum = var_sum + var;
    val_sum = val_sum + value;
    n = n + 1;
%     disp(sqrt(var))
end

val_out = val_sum / n;
var_out = var_sum / n / n;
% disp(['Number of ibslike calls: ' num2str(n)])
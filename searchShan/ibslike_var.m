function [val_out, var_out, n] = ibslike_var(FUN, PARAMS, RESPMAT, DESIGNMAT, options, var_limit)
% function ibslike_var(FUN, PARAMS, RESPMAT, DESIGNMAT, options, var_limit)
% This wrapper function runs ibslike with all parameters just passed
% through until the variance of the estimate is smaller than var_limit

[value, var] = ibslike(FUN, PARAMS, RESPMAT, DESIGNMAT, options);
val_sum = value;
var_sum = var;
n = 1;

while (var_sum / n / n) > var_limit
    [value, var] = ibslike(FUN, PARAMS, RESPMAT, DESIGNMAT, options);
    var_sum = var_sum + var;
    val_sum = val_sum + value;
    n = n + 1;
    % fprintf('%f\n', var_sum / n / n)
end

val_out = val_sum / n;
var_out = var_sum / n / n;
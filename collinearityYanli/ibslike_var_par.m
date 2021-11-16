function [val_out, var_out, n] = ibslike_var_par(FUN, PARAMS, RESPMAT, DESIGNMAT, options, var_limit)
% function ibslike_var(FUN, PARAMS, RESPMAT, DESIGNMAT, options, var_limit)
% This wrapper function runs ibslike with all parameters just passed
% through until the variance of the estimate is smaller than var_limit

[value, var] = ibslike(FUN, PARAMS, RESPMAT, DESIGNMAT, options);
val_sum = value;
var_sum = var;
n = 1;

while (var_sum / n / n) > var_limit
    % go to 90% of the samples we expect to need 
    n_sample = 0.9 *(var_sum / n / var_limit - n);
    % restrict n_sample to the range 1:5n to avoid excessive extrapolation
    n_sample = round(max(min(n_sample, 5*n),1));
    values = zeros(n_sample,1);
    vars = zeros(n_sample,1);
    parfor i = 1:n_sample
        [value, var] = ibslike(FUN, PARAMS, RESPMAT, DESIGNMAT, options);
        values(i) = value;
        vars(i) = var;
    end
    var_sum = var_sum + sum(vars);
    val_sum = val_sum + sum(values);
    n = n + n_sample;
    % fprintf('%f\n', var_sum / n / n)
    % fprintf('%d\n', n_sample);
end

val_out = val_sum / n;
var_out = var_sum / n / n;
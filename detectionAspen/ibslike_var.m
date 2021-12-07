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
%     (var_sum / n / n)
end

val_out = val_sum / n;
var_out = var_sum / n / n;

% BAYES
% ~3:25 pm:  started
%  3:49 pm:  22 fvals
%  3:52 pm:  37 fvals
%  3:56 pm:  66 fvals
%  4:01 pm:  93 fvals
%  4:05 pm: 122 fvals
%  4:10 pm: 152 fvals
%  4:14 pm: 178 fvals
%  4:18 pm: 198 fvals
%  4:23 pm: 239 fvals
%  4:30 pm: 268 fvals
%  4:38 pm: 316
%  4:46 pm: 365
%  5:16 pm: 520
% ~5:48 pm: done!

% FREQ2
%  5:51 pm:  started
%  (took 30 minute nap and closed laptop)
%  6:40 pm:  22 fvals
%  6:43 pm:  33 fvals
%  7:17 pm:  51 fvals
%  7:34 pm:  88 fvals
%  7:40 pm: 109 fvals
%  7:57 pm: 158 fvals
%  8:12 pm: 229 fvals
%  8:44 pm: 365 fvals
%  9:26 pm: 525 fvals



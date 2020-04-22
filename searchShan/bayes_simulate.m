function responses = bayes_simulate(stim, pars)
% function responses = bayes_simulate(stim, pars)
% simulates responses for the Bayes optimal observer
% parameter passed as pars are:
%       sigmaNoise = pars(1)
%       beta = pars(2)
%       bias = pars(3)
%       lapse = pars(4)
%
% Note: In contrast to the original paper we use the log evidence ratio,
% not the difference in probability.

sigmaNoise = pars(1);
beta = pars(2);
bias = pars(3);
lapse = pars(4);

sigma = 9.06;

s = [stim(:,1), stim(:,2), stim(:,2), stim(:,2)];
x = s + sigmaNoise*randn(size(s));

pL = zeros(size(s));

for iL = 1:4
    x_l = x(:, iL);
    x_without_l = [x(:, 1:(iL - 1)),x(:, (iL + 1):end)];
    x_mean_without_l = mean(x_without_l, 2);
    t1 = (x_mean_without_l.^2) ./ 2 ./ (sigma.^2 + sigmaNoise.^2/3);
    t2 = sum((x_without_l-x_mean_without_l).^2 ./ 2 ./ sigmaNoise.^2, 2);
    t3 = x_l.^2 / 2 / (sigmaNoise.^2 + sigma.^2);
    pL(:,iL) = exp(-t1 - t2 - t3);
end
pL = pL ./ sum(pL, 2);
z = sigma .* x ./ sqrt((sigma.^2 + sigmaNoise.^2) .* sigmaNoise.^2);

%p1 = sum(pL .* normcdf(z), 2);
%p2 = sum(pL .* normcdf(-z), 2);
d = log(sum(pL .* normcdf(z), 2)) - log(sum(pL .* normcdf(-z), 2));

decision_variable = bias + beta*d;
p_response = zeros(size(decision_variable));
p_response(decision_variable <=-20) = 0;
p_response(decision_variable >=20) = 1;
selection = (decision_variable>-20) & (decision_variable<20);
p_response(selection) = exp(decision_variable(selection)) ./ (1 + exp(decision_variable(selection)));
p_response = lapse/2 + (1-lapse) * p_response;

responses = rand(size(p_response)) < p_response;


function responses = pe_simulate(stim, pars)
% function responses = pe_simulate(stim, pars)
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

logp = zeros([2,size(s)]);

for iL = 1:4
    x_without_l = [x(:, 1:(iL - 1)),x(:, (iL + 1):end)];
    s_target = (sigma.^2 ./ (sigma.^2 + sigmaNoise.^2)) .* x(:,iL);
    s_target_1 = max(0, s_target);
    s_target_2 = min(0, s_target);
    s_distractor = sum((sigma.^2 ./ (3*sigma.^2 + sigmaNoise.^2)) .* x_without_l, 2);
    logp_dist = -1./2 * (s_distractor.^2 ./ sigma.^2 +  ...
          + sum((x_without_l - s_distractor).^2, 2)./ sigmaNoise.^2);
    logp(2,:,iL) = logp_dist - 1./2 *(s_target_2.^2 ./ sigma.^2 ...
          + (x(:, iL) - s_target_2).^2 ./ sigmaNoise.^2);
    logp(1,:,iL) = logp_dist - 1./2 *(s_target_1.^2 ./ sigma.^2 ...
          + (x(:, iL) - s_target_1).^2 ./ sigmaNoise.^2);
end

% loc = randi(4, size(stim,1), 1);
% [~, loc] = max(reshape(permute(logp, [2, 3, 1]), [], 8), [], 2);
% loc = mod(loc - 1, 4) + 1;
% ind1 = sub2ind(size(logp), ones(size(s, 1), 1), (1:size(s, 1))', loc);
% ind2 = sub2ind(size(logp), 2 * ones(size(s, 1), 1), (1:size(s, 1))', loc);
% d = logp(ind1) - logp(ind2);
d = max(logp(1,:,:),[],3) - max(logp(2,:,:),[],3);
d = d(:);

decision_variable = bias + beta*d;
p_response = zeros(size(decision_variable));
p_response(decision_variable <= -20) = 0;
p_response(decision_variable >= 20) = 1;
selection = (decision_variable > -20) & (decision_variable < 20);
p_response(selection) = exp(decision_variable(selection)) ./ (1 + exp(decision_variable(selection)));
p_response = lapse/2 + (1-lapse) * p_response;

responses = rand(size(p_response)) < p_response;


function responses = sample_simulate(stim, pars)
% function responses = bayes_simulate(stim, pars)
% simulates responses for the Bayes optimal observer
% parameter passed as pars are:
%       sigmaNoise = pars(1)
%       beta = pars(2)
%       bias = pars(3)
%       lapse = pars(4)
%       n_samp = pars(5)
%
% Note: In contrast to the original paper we use the log evidence ratio,
% not the difference in probability.

sigmaNoise = pars(1);
beta = pars(2);
bias = pars(3);
lapse = pars(4);
if length(pars) > 4
    n_samples = floor(pars(5));
    if n_samples > 10000
        n_samples = 10000;
        warning('n_samples reduced to 10000 to avoid extreme computation')
    elseif rand > mod(n_samples, 1)
        n_samples = floor(n_samples);
    else
        n_samples = ceil(n_samples);
    end
else
    n_samples = 1;
end

sigma = 9.06;

s = [stim(:,1), stim(:,2), stim(:,2), stim(:,2)];
x = s + sigmaNoise*randn(size(s));

L1 = randi(4, size(stim,1), n_samples);
L2 = randi(4, size(stim,1), n_samples);
s_samp1 = sigma * randn(size(stim,1), n_samples, 2);
s_samp1(:,:,2) = abs(s_samp1(:,:,2));
s_samp2 = sigma * randn(size(stim,1), n_samples, 2);
s_samp2(:,:,2) = -abs(s_samp2(:,:,2));

idx1 = (L1-1) .* size(stim,1) + (1:size(stim,1))';
idx2 = (L2-1) .* size(stim,1) + (1:size(stim,1))';

% evaluate assumed target position
evidence1 = -(x(idx1) - s_samp1(:,:,2)) .^2 ./ sigmaNoise.^2 ./2;
evidence2 = -(x(idx2) - s_samp2(:,:,2)) .^2 ./ sigmaNoise.^2 ./2;

% evaluate all other positions
for iL = 1:3
    idx1 = mod(L1-1 + iL, 4) .* size(stim, 1) + (1:size(stim, 1))';
    idx2 = mod(L2-1 + iL, 4) .* size(stim, 1) + (1:size(stim, 1))';
    evidence1 = evidence1 -(x(idx1) - s_samp1(:,:,1)) .^2 ./ sigmaNoise.^2 ./2;
    evidence2 = evidence2 -(x(idx2) - s_samp2(:,:,1)) .^2 ./ sigmaNoise.^2 ./2;   
end

d = logsumexp(evidence1, 2) - logsumexp(evidence2, 2);

decision_variable = bias + beta*d;
p_response = zeros(size(decision_variable));
p_response(decision_variable <=-20) = 0;
p_response(decision_variable >=20) = 1;
selection = (decision_variable>-20) & (decision_variable<20);
p_response(selection) = exp(decision_variable(selection)) ./ (1 + exp(decision_variable(selection)));
p_response = lapse/2 + (1-lapse) * p_response;

responses = rand(size(p_response)) < p_response;


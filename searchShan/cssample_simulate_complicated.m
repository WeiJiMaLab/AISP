function responses = cssample_simulate(stim, pars)
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

n_trials = size(stim, 1);

sigma = 9.06;

range = (1:n_trials)';

s = [stim(:,1), stim(:,2), stim(:,2), stim(:,2)];
x = s + sigmaNoise*randn(size(s));

% initialize sampling to prior samples
s1_samp = sigma * randn(n_trials, 2);
s1_samp(:, 1) = abs(s1_samp(:,1));
L1 = randi(4, n_trials, 1);
s2_samp = sigma * randn(n_trials, 2);
s2_samp(:, 1) = -abs(s2_samp(:,1));
L2 = randi(4, n_trials, 1);
% s_1 = nan(n_trials, n_samples, 2);
% s_2 = nan(n_trials, n_samples, 2);
C = nan(n_trials, n_samples);

idx1 = (L1-1) .* n_trials + range;
idx2 = (L2-1) .* n_trials + range;

% evaluate assumed target position
l1_samp = -(x(idx1) - s1_samp(:,1)) .^2 ./ sigmaNoise.^2 ./2;
l2_samp = -(x(idx2) - s2_samp(:,1)) .^2 ./ sigmaNoise.^2 ./2;

% evaluate all other positions
for iL = 1:3
    idx1 = mod(L1-1 + iL, 4) .* n_trials + range;
    idx2 = mod(L2-1 + iL, 4) .* n_trials + range;
    l1_samp = l1_samp -(x(idx1) - s1_samp(:,2)) .^2 ./ sigmaNoise.^2 ./2;
    l2_samp = l2_samp -(x(idx2) - s2_samp(:,2)) .^2 ./ sigmaNoise.^2 ./2;   
end


for i = 1:n_samples
    % sample new s1
    s1_samp_new = sigma * randn(n_trials, 2);
    s1_samp_new(:, 1) = abs(s1_samp_new(:,1));
    L1_new = randi(4, n_trials, 1);
    idx1 = (L1_new-1) .* n_trials + range;
    l1_samp_new = -(x(idx1) - s1_samp_new(:,1)) .^2 ./ sigmaNoise.^2 ./2;
    for iL = 1:3
        idx1 = mod(L1_new-1 + iL, 4) .* n_trials + range;
        l1_samp_new = l1_samp_new -(x(idx1) - s1_samp_new(:,2)) .^2 ./ sigmaNoise.^2 ./2;
    end
    accept = rand(n_trials,1)< exp(l1_samp_new-l1_samp);
    % s1_samp(accept, :) = s1_samp_new(accept, :);
    l1_samp(accept) = l1_samp_new(accept);
    
    % sample s2
    s2_samp_new = sigma * randn(n_trials, 2);
    s2_samp_new(:, 1) = -abs(s2_samp_new(:,1));
    L2_new = randi(4, n_trials, 1);
    idx2 = (L2_new-1) .* n_trials + range;
    l2_samp_new = -(x(idx2) - s2_samp_new(:,1)) .^2 ./ sigmaNoise.^2 ./2;
    for iL = 1:3
        idx2 = mod(L2_new-1 + iL, 4) .* n_trials + range;
        l2_samp_new = l2_samp_new -(x(idx2) - s2_samp_new(:,2)) .^2 ./ sigmaNoise.^2 ./2;   
    end
    accept = rand(n_trials,1)< exp(l2_samp_new-l2_samp);
    % s2_samp(accept, :) = s2_samp_new(accept, :);
    l2_samp(accept) = l2_samp_new(accept);

    % sample C
    p_c = exp(l1_samp-l2_samp) ./ (exp(l1_samp-l2_samp) + 1);
    p_c((l1_samp-l2_samp) > 25) = 1;
    p_c((l1_samp-l2_samp) < -25) = 0;
    C(:,i) = rand(n_trials, 1) < p_c;
    % s_1(:,i,:) = s1_samp;
    % s_2(:,i,:) = s2_samp;
end

d = log(sum(C, 2) + 1) - log(n_samples - sum(C, 2) + 1);
p = lapse/2 + (1-lapse)./(1+exp(bias+beta*d));
responses = (rand(n_trials, 1) < p);
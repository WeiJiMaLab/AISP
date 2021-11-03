function responses = cssample_simulate_complicated(stim, pars)
% function responses = bayes_simulate(data,pars)
% simulates respones of the bayesian observer.
% stim should be a matrix of size trials x 3 with columns:
%     mean 
%     stimulus s1
%     stimulus s2

sigmas = pars(1:4);
beta0 = pars(5);
beta = pars(6);
lambda = pars(7);
if length(pars) > 7
    n_samples = floor(pars(8));
    if n_samples > 10000
        n_samples = 10000;
        warning('n_samples reduced to 10000 to avoid extreme computation')
    end
else
    n_samples = 1;
end
sigma0 = 24;

n_trials = size(stim,1);

s = stim(:, [2,3]) - stim(:, 1);

sigmaNoise = zeros(size(stim, 1), 1);
sigmaNoise(stim(:,1)==0) = sigmas(1);
sigmaNoise(stim(:,1)==240) = sigmas(2);
sigmaNoise(stim(:,1)==480) = sigmas(3);
sigmaNoise(stim(:,1)==840) = sigmas(4);

x = s + repmat(sigmaNoise, [1,2]) .* randn(size(stim, 1), 2);

lsn = log(sigmaNoise(:));
ls0 = log(sigma0);
% initialize sampling to prior samples
s1_samp = sigma0 * repmat(randn(n_trials, 1), [1,2]);
s2_samp = sigma0 * randn(n_trials, 2);
% s_1 = nan(n_trials, n_samples);
% s_2 = nan(n_trials, n_samples);
C = nan(n_trials, n_samples);

l1_samp = - sum((s1_samp-x).^2, 2) ./ 2 ./ sigmaNoise .^2 ...
    - s1_samp(:,1).^2 ./ 2 ./ (sigma0 .^2) - ls0 - lsn;
l2_samp = - sum((s2_samp-x).^2, 2) ./ 2 ./ sigmaNoise .^2 ...
    - sum(s2_samp.^2 ./ 2 ./ (sigma0 .^2), 2) - 2*ls0 - lsn;
% sampling
for i = 1:n_samples
    % sample s1
    s1_samp_new = sigma0 * repmat(randn(n_trials, 1), [1,2]);
    l1_samp_new = - sum((s1_samp-x).^2, 2) ./ 2 ./ sigmaNoise .^2 ...
        - s1_samp(:,1).^2 ./ 2 ./ (sigma0 .^2) - ls0 - lsn;
    accept = rand(n_trials,1)< exp(l1_samp_new-l1_samp);
    % s1_samp(accept) = s1_samp_new(accept);
    l1_samp(accept) = l1_samp_new(accept);
    % sample s2
    s2_samp_new = sigma0 * randn(n_trials, 2);
    l2_samp_new = - sum((s2_samp_new-x).^2, 2) ./ 2 ./ sigmaNoise .^2 ...
        - sum(s2_samp_new.^2 ./ 2 ./ (sigma0 .^2), 2) - 2*ls0 - lsn;
    accept = rand(n_trials, 1)< exp(l2_samp_new-l2_samp);
    % s2_samp(accept) = s2_samp_new(accept);
    l2_samp(accept) = l2_samp_new(accept);
    % sample C
    p_c = exp(l1_samp-l2_samp) ./ (exp(l1_samp-l2_samp) + 1);
    p_c((l1_samp-l2_samp) > 25) = 1;
    p_c((l1_samp-l2_samp) < -25) = 0;
    C(:,i) = rand(n_trials, 1) < p_c;
    % s_1(:,i) = s1_samp;
    % s_2(:,i) = s2_samp;
end
d = log(sum(C, 2) + 1) - log(n_samples - sum(C, 2) + 1);
p = lambda/2 + (1-lambda)./(1+exp(beta0+beta*d));
responses = (rand(n_trials, 1) < p);

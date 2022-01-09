function responses = cssample_simulate(stim, pars)
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
    n_samples = pars(8);
    % random averaging between n_samples values
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
sigma0 = 24;

n_trials = size(stim,1);

s = stim(:, [2,3]) - stim(:, 1);

sigmaNoise = zeros(size(stim, 1), 1);
sigmaNoise(stim(:,1)==0) = sigmas(1);
sigmaNoise(stim(:,1)==240) = sigmas(2);
sigmaNoise(stim(:,1)==480) = sigmas(3);
sigmaNoise(stim(:,1)==840) = sigmas(4);

x = s + repmat(sigmaNoise, [1,2]) .* randn(size(stim, 1), 2);

% initialize sampling to prior samples
s_samp = sigma0 * randn(n_trials, 2);
C_samp = randi(2, n_trials, 1);
s_samp(C_samp==1, 2) = s_samp(C_samp==1,1);
% s = nan(n_trials, n_samples);
C = nan(n_trials, n_samples);


l_samp = - sum((s_samp-x).^2, 2) ./ 2 ./ sigmaNoise .^2;
% sampling
for i = 1:n_samples
    % sample s1
    s_samp_new = sigma0 * randn(n_trials, 2);
    C_samp_new = randi(2, n_trials, 1);
    s_samp_new(C_samp_new==1, 2) = s_samp_new(C_samp_new==1,1);
    l_samp_new = - sum((s_samp_new-x).^2, 2) ./ 2 ./ sigmaNoise .^2;
    accept = rand(n_trials,1)< exp(l_samp_new-l_samp);
    % s_samp(accept) = s_samp_new(accept);
    C_samp(accept) = C_samp_new(accept);
    l_samp(accept) = l_samp_new(accept);
    C(:,i) = C_samp;
    % s(:,i) = s_samp;
end
d = log(sum(C==1, 2) + 1) - log(n_samples - sum(C==1, 2) + 1);
p = lambda/2 + (1-lambda)./(1+exp(-beta0-beta*d));
responses = (rand(n_trials, 1) < p);

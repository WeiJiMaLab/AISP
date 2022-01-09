function responses = cssample_simulate(data,sigmas,beta,lambda,n_samples)
% function responses = samp_simulate(data,sigmas,beta,lambda)
% simulates respones of the second point estimate observer.
% data should be a matrix of size trials x 2 with columns:
%     reliability [1-6] 
%     stimulus s 
% sigmas should be a 6 element vector of noise standard deviations
% beta should be the constant and slope for the logistic link
% lambda is the lapserate

if ~exist('lambda','var') || isempty(lambda)
    lambda = 0;
end
if ~exist('n_samples','var') || isempty(n_samples)
    n_samples = 1;
else
    n_samples = floor(n_samples);
end
if n_samples > 10000
    n_samples = 10000;
    warning('n_samples reduced to 10000 to avoid extreme computation')
end

sigma1 = 3;
sigma2 = 12;
sigma12 = [sigma1, sigma2];
n_trials = size(data,1);
sigmaNoise = sigmas(floor(data(:,1)))';

s = data(:,2);
x = s + sigmaNoise(:) .* randn(n_trials, 1);
lsn = log(sigmaNoise(:));
% initialize sampling to prior samples
C_samp = randi(2, n_trials, 1);
s_samp = sigma12(C_samp)' .* randn(n_trials, 1);
s_s = nan(n_trials, n_samples);
C = nan(n_trials, n_samples);

l_samp = - (s_samp-x).^2 ./ 2 ./ sigmaNoise .^2 - lsn;
% sampling
for i = 1:n_samples
    % sample C
    C_new = randi(2, n_trials, 1);
    % sample s
    s_samp_new = sigma12(C_new)' .* randn(n_trials, 1);
    l_samp_new = - (s_samp_new-x).^2 ./ 2 ./ sigmaNoise .^2 - lsn;
    accept = rand(n_trials,1) < exp(l_samp_new-l_samp);
    s_samp(accept) = s_samp_new(accept);
    l_samp(accept) = l_samp_new(accept);
    C_samp(accept) = C_new(accept);
    C(:,i) = C_samp;
    s_s(:,i) = s_samp;
end
% 1 is the wide category = 2 here
d = log(sum(C==1, 2) + 1) - log(sum(C==2, 2) + 1);
p = lambda/2 + (1-lambda)./(1+exp(beta(1)+beta(2)*d));
responses = (rand(n_trials, 1) < p);
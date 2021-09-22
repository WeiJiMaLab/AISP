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
n_trials = size(data,1);
sigmaNoise = sigmas(floor(data(:,1)))';

s = data(:,2);
x = s + sigmaNoise(:) .* randn(n_trials, 1);
lsn = log(sigmaNoise(:));
ls1 = log(sigma1);
ls2 = log(sigma2);
% initialize sampling to prior samples
s1_samp = sigma1 * randn(n_trials, 1);
s2_samp = sigma2 * randn(n_trials, 1);
% s_1 = nan(n_trials, n_samples);
% s_2 = nan(n_trials, n_samples);
C = nan(n_trials, n_samples);

l1_samp = - (s1_samp-x).^2 ./ 2 ./ sigmaNoise .^2 ...
    - s1_samp.^2 ./ 2 ./ (sigma1 .^2) - ls1 - lsn;
l2_samp = - (s2_samp-x).^2 ./ 2 ./ sigmaNoise .^2 ...
    - s2_samp.^2 ./ 2 ./ (sigma2 .^2) - ls2 - lsn;
% sampling
for i = 1:n_samples
    % sample s1
    s1_samp_new = sigma1 * randn(n_trials,1);
    l1_samp_new = - (s1_samp_new-x).^2 ./ 2 ./ sigmaNoise .^2 ...
        - s1_samp_new.^2 ./ 2 ./ (sigma1 .^2) - ls1 - lsn;
    accept = rand(n_trials,1)< exp(l1_samp_new-l1_samp);
    % s1_samp(accept) = s1_samp_new(accept);
    l1_samp(accept) = l1_samp_new(accept);
    % sample s2
    s2_samp_new = sigma2 * randn(n_trials,1);
    l2_samp_new = - (s2_samp_new-x).^2 ./ 2 ./ sigmaNoise .^2 ...
        - s2_samp_new.^2 ./ 2 ./ (sigma2 .^2) - ls2 - lsn;
    accept = rand(n_trials,1)< exp(l2_samp_new-l2_samp);
    % s2_samp(accept) = s2_samp_new(accept);
    l2_samp(accept) = l2_samp_new(accept);
    % sample C
    p_c = exp(l1_samp-l2_samp) ./ (exp(l1_samp-l2_samp) + 1);
    p_c((l1_samp-l2_samp)>25) = 1;
    p_c((l1_samp-l2_samp)<-25) = 0;
    C(:,i) = rand(n_trials, 1) < p_c;
    % s_1(:,i) = s1_samp;
    % s_2(:,i) = s2_samp;
end
% 1 is the wide category = 2 here
d = log(n_samples - sum(C, 2) + 1) - log(sum(C, 2) + 1);
p = lambda/2 + (1-lambda)./(1+exp(beta(1)+beta(2)*d));
responses = (rand(n_trials, 1) < p);
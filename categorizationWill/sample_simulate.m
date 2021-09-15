function responses = sample_simulate(data,sigmas,beta,lambda,n_samples)
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
sigmaNoise = sigmas(floor(data(:,1)));
s = data(:,2);
x = s + sigmaNoise(:) .* randn(n_trials, 1);
s1_samp = sigma1 * randn(n_trials, n_samples);
s2_samp = sigma2 * randn(n_trials, n_samples);
evidence1 = logsumexp(- (s1_samp-repmat(x, 1, n_samples)).^2 ./ 2 ./ repmat(sigmaNoise' .^2, 1, n_samples), 2);
evidence2 = logsumexp(- (s2_samp-repmat(x, 1, n_samples)).^2 ./ 2 ./ repmat(sigmaNoise' .^2, 1, n_samples), 2);
d = evidence1 - evidence2;
p = lambda/2 + (1-lambda)./(1+exp(beta(1)+beta(2)*d));
responses = (rand(n_trials, 1) < p);
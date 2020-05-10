function responses = pe_simulate(stim, pars)
% function responses = bayes_simulate(data,sigmas,beta,lambda)
% simulates respones of the bayesian observer.
% stim should be a matrix of size trials x 3 with columns:
%     mean 
%     stimulus s1
%     stimulus s2 
% sigmas should be a 3 element vector of noise standard deviations
% beta should be the constant and slope for the logistic link
% lambda is the lapserate


sigmas = pars(1:4);
beta0 = pars(5);
beta = pars(6);
lambda = pars(7);
sigma0 = 24;

s = stim(:, [2,3]) - stim(:, 1);

sigmaNoise = zeros(size(stim, 1), 1);
sigmaNoise(stim(:,1)==0) = sigmas(1);
sigmaNoise(stim(:,1)==240) = sigmas(2);
sigmaNoise(stim(:,1)==480) = sigmas(3);
sigmaNoise(stim(:,1)==840) = sigmas(4);

x = s + repmat(sigmaNoise, [1,2]) .* randn(size(stim, 1), 2);

shatSeparate = (sigma0.^2 ./ (sigma0.^2 + sigmaNoise.^2)) .* x;
shatSame = (sigma0.^2 ./ (2*sigma0.^2 + sigmaNoise.^2)) .* sum(x, 2);

d = - 1 / 2 .* ((sum(shatSeparate.^2, 2) - shatSame.^2) ./ sigma0.^2 ...
    + (sum((x - shatSeparate).^2, 2) - sum((x - repmat(shatSame,[1,2])).^2, 2))./ sigmaNoise);
p = lambda./2 + (1-lambda)./(1+exp(beta0 + beta*d));
responses = (rand(size(stim, 1), 1) < p);

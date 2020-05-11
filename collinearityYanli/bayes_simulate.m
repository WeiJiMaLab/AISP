function responses = bayes_simulate(stim, pars)
% function responses = bayes_simulate(data,pars, exact)
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

priorD = sigmas .* sqrt(sigmas.^2 + 2*sigma0.^2) ./ (sigmas.^2 + sigma0.^2);
priorD = log(priorD);
        
s = stim(:, [2,3]) - stim(:, 1);

sigmaNoise = zeros(size(stim, 1), 1);
sigmaNoise(stim(:,1)==0) = sigmas(1);
sigmaNoise(stim(:,1)==240) = sigmas(2);
sigmaNoise(stim(:,1)==480) = sigmas(3);
sigmaNoise(stim(:,1)==840) = sigmas(4);
pd = zeros(size(stim, 1), 1);
pd(stim(:,1)==0) = priorD(1);
pd(stim(:,1)==240) = priorD(2);
pd(stim(:,1)==480) = priorD(3);
pd(stim(:,1)==840) = priorD(4);

x = s + repmat(sigmaNoise, [1,2]) .* randn(size(stim, 1), 2);
% old wrong formulas:
%d = pd - sigma0.^2 ./ 2 ./ sigmaNoise.^2 ./ (2 .* sigmaNoise.^2 + sigma0.^2) .* x(:,1) .* x(:,2);
%d = d + sigma0.^2 ./ (sigma0.^2 + sigmaNoise.^2) ./ (4 .* sigmaNoise.^2 + 2 * sigma0.^2) .* (x(:,1).^2 + x(:,2).^2);
d = pd - (x(:,1).^2 + x(:,2).^2)./2./(sigma0.^2 + sigmaNoise.^2);
d = d  - (x(:,1) + x(:,2)).^2 .* sigma0.^2 ./ 2./ sigmaNoise.^2 ./ (2*sigma0.^2 + sigmaNoise.^2);
d = d  + (x(:,1).^2 + x(:,2).^2)./2 ./ sigmaNoise.^2;
p = lambda./ 2 + (1-lambda)./(1+exp(beta0 + beta.*d));
responses = (rand(size(stim, 1), 1) < p);

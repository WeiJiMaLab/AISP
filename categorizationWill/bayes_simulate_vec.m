function responses = bayes_simulate_vec(data,sigmas,beta,lambda)
% function responses = bayes_simulate(data,sigmas,beta,lambda)
% simulates respones of the bayesian observer.
% data should be a matrix of size trials x 2 with columns:
%     reliability [1-6] 
%     stimulus s 
% sigmas should be a 6 element vector of noise standard deviations
% beta should be the constant and slope for the logistic link
% lambda is the lapserate

if ~exist('lambda','var') || isempty(lambda)
    lambda = 0;
end
   
sigma1 = 3;
sigma2 = 12;
sigmaNoise = sigmas(floor(data(:,1)))';
s = data(:,2);
x = s + sigmaNoise.*randn(length(s),1);
d = 1/2*log((sigma2.^2+sigmaNoise.^2)./(sigma1.^2+sigmaNoise.^2))- ...
    x.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2+sigmaNoise.^2)./(sigma2^2+sigmaNoise.^2);
p = lambda/2 + (1-lambda)./(1+exp(beta(1)+beta(2).*d));
responses = (rand(length(p),1)<p);
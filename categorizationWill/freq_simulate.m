function responses = freq_simulate(data,sigmas,beta,lambda)
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
responses = zeros(size(data,1),1);
for iTrial = 1:size(data,1)
    sigmaNoise = sigmas(floor(data(iTrial,1)));
    s = data(iTrial,2);
    x = s + sigmaNoise*randn;
    d = 1/2*log((sigma2.^2)./(sigma1.^2))- ...
        x.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2+sigmaNoise.^2)./(sigma2^2+sigmaNoise^2);
    p = lambda/2 + (1-lambda)/(1+exp(beta(1)+beta(2)*d));
    responses(iTrial) = (rand<p);
end
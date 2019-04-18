function responses = freq3_simulate(data,sigmas,beta,lambda)
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
    w1 = normpdf(x,0,sigma1^2+sigmaNoise);
    w2 = normpdf(x,0,sigma2^2+sigmaNoise);
    w1n = w1./(w1+w2);
    w2n = w2./(w1+w2);
    shat = w1n .* sigma1./(sigma1+sigmaNoise).*x+w2n .* sigma2./(sigma2+sigmaNoise).*x;
    d = 1/2*log((sigma2.^2)./(sigma1.^2))- ...
        shat.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2)./(sigma2^2);
    p = lambda/2 + (1-lambda)/(1+exp(beta(1)+beta(2)*d));
    responses(iTrial) = (rand<p);
end
function likelihood = freq3_likelihood(data,sigmas,beta,lambda,Nsamples)
% function likelihood = bayes_likelihood(data,sigmas,beta,lambda)
% computes the likelihood for the bayesian observer.
% data should be a matrix of size trials x 3 with columns:
%     reliability [1-6] 
%     stimulus s 
%     reponse [1,-1]
% sigmas should be a 6 element vector of noise standard deviations
% beta should be the constant and slope for the logistic link
% lambda is the lapserate

if ~exist('lambda','var') || isempty(lambda)
    lambda = 0;
end
   
if ~exist('Nsamples','var') || isempty(Nsamples)
    Nsamples = 1000;
end
sigma1 = 3;
sigma2 = 12;
likelihood = zeros(size(data,1),1);
%ds = zeros(size(data,1),1);
%ps = zeros(size(data,1),1);
for iTrial = 1:size(data,1)
    sigmaNoise = sigmas(data(iTrial,1));
    s = data(iTrial,2);
    x = s + sigmaNoise*randn(Nsamples,1);
    w1 = normpdf(x,0,sigma1^2+sigmaNoise);
    w2 = normpdf(x,0,sigma2^2+sigmaNoise);
    w1n = w1./(w1+w2);
    w2n = w2./(w1+w2);
    shat = w1n .* sigma1./(sigma1+sigmaNoise).*x+w2n .* sigma2./(sigma2+sigmaNoise).*x;
    d = 1/2*log((sigma2.^2)./(sigma1.^2))- ...
        shat.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2)./(sigma2^2);
    p = mean(lambda/2 + (1-lambda)./(1+exp(-beta(1)-beta(2).*d)));
    %ds(iTrial) = mean(d);
    %ps(iTrial) = p;
    likelihood(iTrial) = (data(iTrial,3)==-1)*p + (data(iTrial,3)==1)*(1-p);
end
function responses = var_simulate(data,sigmas,beta,lambda)
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
sigma12 = sigma1.^2;
sigma22 = sigma2.^2;
sigm = sigmas(data(:,1));
tol = 10^-5;
responses = zeros(size(data,1),1);
for iTrial = 1:size(data,1)
    sigmaNoise = sigmas(floor(data(iTrial,1)));
    s = data(iTrial,2);
    x = s + sigmaNoise*randn;
    q1 = zeros(size(x));
    q1new = 0.5 * ones(size(x));
    k = 0;
    while any(abs(q1new-q1)>tol) && ~(k>50)
        k = k+1;
        q1 = q1new;
        sigmanew = (sigma12.*sigma22)./(q1*sigma12.*sigmaNoise.^2 + ...
            (1-q1).*sigma2.^2.*sigmaNoise.^2+sigma1.^2*sigma2.^2);
        shat = sigmanew.*x;
        svar = sigmanew.*sigmaNoise.^2;
        z = -(shat.^2+svar)./2;
        z1 = sigma2.*exp(z./sigma12);
        z2 = sigma1.*exp(z./sigma22);
        q1new = z1./(z1 + z2);
    end
    q1 = q1new;
    d = log(q1)-log(1-q1);
    p = lambda/2 + (1-lambda)/(1+exp(beta(1)+beta(2)*d));
    responses(iTrial) = (rand<p);
end
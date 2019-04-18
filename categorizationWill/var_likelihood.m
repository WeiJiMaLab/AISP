function likelihood = var_likelihood(data,sigmas,beta,lambda,Nsamples)
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
tol = 0.00001;
likelihood = zeros(size(data,1),1);

%ds = zeros(size(data,1),1);
%ps = zeros(size(data,1),1);
sigm = sigmas(data(:,1));
sigma12 = sigma1.^2;
sigma22 = sigma2.^2;

for iTrial = 1:size(data,1)
    sigmaNoise = sigm(iTrial);
    s = data(iTrial,2);
    x = s + sigmaNoise*randn(Nsamples,1);
    q1 = zeros(size(x));
    q1new = 0.5 * ones(size(x));
    while any(abs(q1new-q1)>tol)
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
    p = mean(lambda/2 + (1-lambda)./(1+exp(-beta(1)-beta(2).*d)));
    %ds(iTrial) = mean(d);
    %ps(iTrial) = p;
    likelihood(iTrial) = (data(iTrial,3)==-1)*p + (data(iTrial,3)==1)*(1-p);
end

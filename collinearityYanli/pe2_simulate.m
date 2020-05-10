function responses = pe2_simulate(stim, pars)
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

responses = zeros(size(stim, 1), 1);
for iTrial = 1:size(stim,1)
    if stim(iTrial,1) == 0
        sigmaNoise = sigmas(1);
    elseif stim(iTrial,1) == 240
        sigmaNoise = sigmas(2);
    elseif stim(iTrial,1) == 480
        sigmaNoise = sigmas(3);
    elseif stim(iTrial,1) == 840
        sigmaNoise = sigmas(4);
    end
    s = stim(iTrial, [2,3]) - stim(iTrial, 1);
    x = s + sigmaNoise * randn(1,2);
    shatSeparate = (sigma0.^2 / (sigma0.^2 + sigmaNoise.^2)) .* x;
    shatSame = (sigma0.^2 / (2*sigma0.^2 + sigmaNoise.^2)) .* sum(x);
    d = - 1 / 2 * ((sum(shatSeparate.^2) - shatSame.^2) ./ sigma0.^2 ...
        + (sum((x - shatSeparate).^2) - sum((x - shatSame).^2))./ sigmaNoise);
    p = lambda/2 + (1-lambda)/(1+exp(beta0 + beta*d));
    responses(iTrial) = (rand < p);
end
function lik = likelihood_optim(data,pars,type)
% function likelihood = bayes_likelihood(data,sigmas,beta,lambda)
% computes the likelihood for the bayesian observer.
% data should be a matrix of size trials x 3 with columns:
%     reliability [1-6] 
%     stimulus s 
%     reponse [1,-1]
% pars should be a 9 element vector

sigmas = exp(pars(1:6));
beta = pars(7:8);
beta(2) = exp(beta(2));
lambda = pars(9);
switch type
    case 'bayes'
        likelihood = bayes_likelihood(data(:,2:4),sigmas,beta,lambda);
    case 'freq'
        likelihood = freq_likelihood(data(:,2:4),sigmas,beta,lambda);
    case 'freq2'
        likelihood = freq2_likelihood(data(:,2:4),sigmas,beta,lambda);
    case 'freq3'
        likelihood = freq3_likelihood(data(:,2:4),sigmas,beta,lambda);
    case 'var'
        likelihood = var_likelihood(data(:,2:4),sigmas,beta,lambda);
end
lik = -sum(log(likelihood));
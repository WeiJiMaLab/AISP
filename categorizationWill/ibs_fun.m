function response = ibs_fun(data,pars,type)
% function response = ibs_fun(data,pars,type)
% simulates data for ibs-sampling
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
        response = bayes_simulate(data(:,2:4),sigmas,beta,lambda);
    case 'freq'
        response = freq_simulate(data(:,2:4),sigmas,beta,lambda);
    case 'freq2'
        response = freq2_simulate(data(:,2:4),sigmas,beta,lambda);
    case 'freq3'
        response = freq3_simulate(data(:,2:4),sigmas,beta,lambda);
    case 'var'
        response = var_simulate(data(:,2:4),sigmas,beta,lambda);
    case 'sample'
        n_samples = pars(10);
        response = sample_simulate(data(:,2:4),sigmas,beta,lambda, n_samples);
    case 'sample1'
        response = sample_simulate(data(:,2:4),sigmas,beta,lambda, 1);
end
response = -1+2*response; % convert to -1,1
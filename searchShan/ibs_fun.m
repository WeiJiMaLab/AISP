function response = ibs_fun(data,pars,type)
% function response = ibs_fun(data,pars,type)
% simulates data for ibs-sampling
% data should be a matrix of size trials x (2 or 3) with columns:
%     stimulus s target
%     stimulus s distractor
%     reponse [1,-1] [optional]
% pars should be a 4 element vector [sigmaNoise, beta, bias, lapse]

switch type
    case 'bayes'
        response = bayes_simulate(data, pars);
    case 'freq'
        response = pe_simulate(data, pars);
end
response = -1+2*response; % convert to -1,1
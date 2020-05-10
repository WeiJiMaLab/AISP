function response = ibs_fun(data, pars, type)
% function response = ibs_fun(data,pars,type)
% simulates data for ibs-sampling
% data can 

switch type
    case 'bayes'
        response = bayes_simulate(data(:,[2,5,6]), pars);
    case 'freq'
        response = pe_simulate(data(:,[2,5,6]), pars);
    case 'freq2'
        response = pe2_simulate(data(:,[2,5,6]), pars);
end
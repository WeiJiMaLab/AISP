function response = ibs_fun(data, pars, type)
% function response = ibs_fun(data,pars,type)
% simulates data for ibs-sampling
% data can 

switch type
    case 'bayes'
        response = bayes_simulate(data(:,[2,5,6]), pars);
    case 'PE'
        response = pe_simulate(data(:,[2,5,6]), pars);
    case 'PE2'
        response = pe2_simulate(data(:,[2,5,6]), pars);
    case {'sample', 'sample1'}
        response = sample_simulate(data(:,[2,5,6]), pars);
end
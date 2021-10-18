function Config = produceConfig()

% 5 models
Config.ModelLabel = {'Bayes', 'Point estimate (Dep. L)', ...
    'Opt. point estimate', 'Point estimate (Ind. L)', ...
    'Importance sampling', 'Joint posterior sampling'};
Config.ModelList =  {'Bayes', 'PE', 'PE2', 'PE_imagineL', ...
    'impSamp', 'jointPostSamp'};
Config.Nreps = 20;

% 3 models
% Config.ModelLabel = {'Bayes', 'Point estimate', 'Opt. point estimate'};
% Config.ModelList =  {'Bayes', 'PE', 'PE2'};
% Config.Nreps = 20;

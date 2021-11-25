function Config = produceConfig()

% 6 models
Config.ModelLabel = {'Bayes', 'Point estimate (Dep. L)', ...
    'Opt. point estimate', 'Point estimate (Ind. L)', ...
    'Importance sampling', 'Joint posterior sampling'};
Config.ModelList =  {'Bayes', 'PE', 'PE2', 'PE_imagineL', ...
    'impSamp', 'jointPostSamp'};
Config.Nreps = 1;
Config.FracPtpnts = 3/13; % Fraction of participants to run

% 3 models
% Config.ModelLabel = {'Bayes', 'Point estimate', 'Opt. point estimate'};
% Config.ModelList =  {'Bayes', 'PE', 'PE2'};
% Config.Nreps = 20;

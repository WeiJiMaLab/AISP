function Config = produceConfig(varargin)

% varargin{1}: bool. If true, don't plot the additional point estimate 
%   observer model. Default false.

if (length(varargin) > 0) && (~isempty(varargin{1}))
    exclExtraPE = varargin{1};
else
    exclExtraPE = false;
end

if ~exclExtraPE
    % 6 models
    Config.ModelLabel = {'Bayes', 'Point estimate (Dep. L)', ...
        'Opt. point estimate', 'Point estimate (Ind. L)', ...
        'Importance sampling', 'Joint posterior sampling'};
    Config.ModelList =  {'Bayes', 'PE', 'PE2', 'PE_imagineL', ...
        'impSamp', 'jointPostSamp'};
    Config.Nreps = 20;
    Config.FracPtpnts = 1; % Fraction of participants to run

elseif exclExtraPE
    % 5 models
    Config.ModelLabel = {'Bayes', 'Point estimate', ...
        'Opt. point estimate', ...
        'Importance sampling', 'Joint posterior sampling'};
    Config.ModelList =  {'Bayes', 'PE', 'PE2', ...
        'impSamp', 'jointPostSamp'};
    Config.Nreps = 20;
    Config.FracPtpnts = 1; % Fraction of participants to run

end

% 3 models
% Config.ModelLabel = {'Bayes', 'Point estimate', 'Opt. point estimate'};
% Config.ModelList =  {'Bayes', 'PE', 'PE2'};
% Config.Nreps = 20;

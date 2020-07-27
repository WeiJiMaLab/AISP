function bestPars = aisp_collectBestFittingParams(nLogLs, pars)
% Finds the parameters from the fit repetition that resulted in the greatest
% loglikelihood

% INPUT
% nLogLs: [num participants x num repetitions] array describing the negative
% log-likelihood resulting from different fitting repititions for a single model
% pars:  [num participants x num params x num repetitions] array

% OUTPUT
% bestPars: [num participants x num params] array

[~, bestRep] = min(nLogLs, [], 2);

Nptpnts = size(pars, 1);
Nparams = size(pars, 2);
bestPars = nan(Nptpnts, Nparams);

for iP = 1 : Nptpnts
    bestPars(iP, :) = pars(iP, :, bestRep(iP));
end



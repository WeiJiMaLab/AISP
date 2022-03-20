function compareParamsAndGoodnessOfFit(DSet, goodnessModel, paramsModel)
% Plot paramter fits for one model against maximum LL of another (or the
% same) model

% INPUT
% DSet: Fitted dataset including parameters and LLs, in the standard format
% goodnessModel: The number of the model, as ordered in DSet, to use for
%   the log-likeliohood (LL) values
% paramsModel: The number of the model, as ordered in DSet, to use for
%   the fitted parameters

% NOTES
% Does not plot the sampling parameters

% JCT, 2022

models = mT_findAppliedModels(DSet);

figure; subplot(4, 1, 1)
hold on
for iP = 1 : length(DSet.P(1).Models(paramsModel).BestFit.Params.LnKappa_x)
    scatter(mT_stackData(...
        DSet.P, @(St)St.Models(goodnessModel).BestFit.LL), ...
        mT_stackData(DSet.P, ...
            @(St)St.Models(paramsModel).BestFit.Params.LnKappa_x(iP)))
end
title('LnKappa_x')

subplot(4, 1, 2)
scatter(mT_stackData(...
    DSet.P, @(St)St.Models(goodnessModel).BestFit.LL), ...
    mT_stackData(DSet.P, ...
    @(St)St.Models(paramsModel).BestFit.Params.Beta0))
title('Beta0')

subplot(4, 1, 3)
scatter(mT_stackData(...
    DSet.P, @(St)St.Models(goodnessModel).BestFit.LL), ...
    mT_stackData(DSet.P, ...
    @(St)St.Models(paramsModel).BestFit.Params.Beta1))
title('Beta1')

subplot(4, 1, 4)
scatter(mT_stackData(...
    DSet.P, @(St)St.Models(goodnessModel).BestFit.LL), ...
    mT_stackData(DSet.P, ...
    @(St)St.Models(paramsModel).BestFit.Params.LapseRate))
title('LapseRate')

for iP = 1 : 4
    subplot(4, 1, iP)
    xlabel(['Log-likelihood; Model: ', models{goodnessModel}])
    ylabel({'Param value', ['Model: ', models{paramsModel}]})
end
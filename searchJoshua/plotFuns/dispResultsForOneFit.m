function dispResultsForOneFit(DSet, parsDir, participant, model, configFile)
% Compare the results of different runs for the same fit

% INPUT
% DSet: Produced by function 'convertToDSetFormat'
% parsDir: Directory containing the fitted parameter data from the cluster
% participant: Number of participant to look at
% model: Number of model to look at
% configFile: string. File path to matlab file to be loaded.

% Collect some info
[X0,LB,UB,PLB,PUB] = get_bads_bounds();
numParams = length(X0);
numPlots = numParams + 1;

Config = load(configFile);
savedName = [parsDir '\pars_' Config.ModelList{model}];
LoadedData = load([savedName,'.mat']);
fittedPars = LoadedData.pars;
fittedPars = fittedPars(participant, :, :);

currentPlot = 1;

%% Likelihood plot

allLLs = mT_stackData(DSet.P(participant).Models(model).Fits, @(st)st.LL);
[maxLL, maxIndex] = max(allLLs);

subplot(ceil(numPlots/2), 2, currentPlot)
histogram(allLLs)
xline(maxLL, 'LineWidth', 3)
xlabel('LL')

currentPlot = currentPlot + 1;

%% Plots for the parameters

assert(isequal(DSet.P(participant).Models(model).Settings.ModelName, ...
    Config.ModelList{model}))

paramNames = {'LnKappa_x', 'LnKappa_x', 'LnKappa_x', 'LnKappa_x', 'Beta0', 'Beta1', 'LapseRate'};
warning('Param names hard coded in. Change this.')

for iPm = 1 : numParams
    paramVals(:) = fittedPars(1, iPm, :);
    
    subplot(ceil(numPlots/2), 2, currentPlot)
    histogram(paramVals)
    
    xline(paramVals(maxIndex), 'LineWidth', 3)
    xlabel(paramNames{iPm})
    
    currentPlot = currentPlot + 1;
end

    
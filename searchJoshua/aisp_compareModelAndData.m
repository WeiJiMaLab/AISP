function Figures = aisp_compareModelAndData(dataDir, parsDir, model)

[SimDSet, DSet] = aisp_simDataset(dataDir, parsDir, model);

% Add summary stats to the dataset
% TODO remove the duplication here
for iPtpnt = 1 : length(SimDSet.P)  
    SimDSet.P(iPtpnt).Data ...
        = computeStimStats(SimDSet.P(iPtpnt).Data, 'circ', false);
end
for iPtpnt = 1 : length(DSet.P)  
    DSet.P(iPtpnt).Data ...
        = computeStimStats(DSet.P(iPtpnt).Data, 'circ', false);
end

Settings.ModelPlotType = 'errorShading';
Settings.DataPlotType = 'scatter';
Settings.NumBins = 10;

Figures = plotCompareTwoDatasets(DSet, SimDSet, Settings);
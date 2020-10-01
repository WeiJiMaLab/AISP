function [bigPlot, indiviudalPlots] = plotAllModelFits(dataDir, parsDir)
% Makes plots using plotHitAndFaAllStats, one for each model, comparing the
% model and data, and then combines them all into one big plot. Also returns the
% individual plots.

Config = load('config.mat');


%% Make the individual plots
smallPlots = cell(length(Config.ModelList), 1);

for iM = 1 : length(Config.ModelList)
    [DSet, SimDSet] = aisp_simDataset(dataDir, parsDir, Config.ModelList{iM});
    
    smallPlots{iM} = figure;
    smallPlots{iM} = plotHitAndFaAllStats(DSet, smallPlots{iM}, 'scatter');
    smallPlots{iM} = plotHitAndFaAllStats(SimDSet, smallPlots{iM}, 'errorShading');
end

indiviudalPlots = smallPlots;


%% Combine the plots

numPlotsInEachSmall = 3;
bigPlot = figure;

for iM = 1 : length(Config.ModelList)
    for iPlotCol = 1 : numPlotsInEachSmall
        bigPlotAx{iM, iPlotCol} = subplot(length(Config.ModelList), ...
            numPlotsInEachSmall, ...
            ((numPlotsInEachSmall*(iM-1)) + iPlotCol));
    end
end

for iM = 1 : length(Config.ModelList)
    smallPlotAllAx = findall(smallPlots{iM}, 'type', 'axes');
    assert(length(smallPlotAllAx) == numPlotsInEachSmall)
    
    % Copy the contents of each small plot
    copied = cell(numPlotsInEachSmall);
    for iPlotCol = 1 : numPlotsInEachSmall
        
        copied{iPlotCol} = copyobj( ...
            smallPlotAllAx(1+numPlotsInEachSmall-iPlotCol), bigPlot);
        copied{iPlotCol}.Position = bigPlotAx{iM, iPlotCol}.Position;
        
        % Some properties need changing for the big plot
        if iM ~= length(Config.ModelList)
            xlabel(copied{iPlotCol}, '')
        end
        
        if iPlotCol == 1
            if iM == 2
               firstLabel = [currentLabel{:}];
            else
                firstLabel = ' ';
            end
            
            currentLabel = copied{iPlotCol}.YLabel.String;
            ylabel(copied{iPlotCol}, ...
                {['{\bf ' firstLabel '}'], Config.ModelLabel{iM}});
        end
    end
end


for iAx = 1 : length(bigPlotAx(:))
    delete(bigPlotAx{iAx})
end
        

% TODO
warning(['Code does not yet remove unecessary numbers for axes, where ', ...
    'axes are duplicated on several rows/columns.'])
        
        
        
        
        

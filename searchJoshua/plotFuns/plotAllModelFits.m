function [bigPlot, indiviudalPlots] = plotAllModelFits(dataDir, parsDir, ...
    configFile)
% Makes plots using plotHitAndFaAllStats, one for each model, comparing the
% model and data, and then combines them all into one big plot. Also returns the
% individual plots.

% INPUT
% configFile: string. File path to matlab file to be loaded.

Config = load(configFile);


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

% Each row/column will have a different offset 
% TODO this is hard coded -- the number of rows and columns. Change
xOffset = [0, -0.03, -0.06, -0.09];
yOffset = [0, 0.03, 0.06, 0.09];

for iM = 1 : length(Config.ModelList)
    for iPlotCol = 1 : numPlotsInEachSmall
        bigPlotAx{iM, iPlotCol} = subplot(length(Config.ModelList), ...
            numPlotsInEachSmall, ...
            ((numPlotsInEachSmall*(iM-1)) + iPlotCol));
        
        set(bigPlotAx{iM, iPlotCol}, 'Units', 'normalized');
        pos = get(bigPlotAx{iM, iPlotCol}, 'Position');
        pos(1) = pos(1) + xOffset(iPlotCol);
        pos(2) = pos(2) + yOffset(iM);
        set(bigPlotAx{iM, iPlotCol}, 'Position', pos);
    end
end

for iM = 1 : length(Config.ModelList)
    smallPlotAllAx = findall(smallPlots{iM}, 'type', 'axes');
    assert(length(smallPlotAllAx) == numPlotsInEachSmall)
    
    % Copy the contents of each small plot
    copied = cell(numPlotsInEachSmall, 1);
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
        
        % Removing duplicated tick labels
        if iM < length(Config.ModelList)
            set(copied{iPlotCol}, 'XTickLabel', [])
        end
        
        if iPlotCol > 1
            set(copied{iPlotCol}, 'YTickLabel', [])
        end
    end
end


for iAx = 1 : length(bigPlotAx(:))
    delete(bigPlotAx{iAx})
end
        
        
        
        
        

function plotFig = aisp_plotOverallPerformance(DSet, varargin)

% INPUT
% varargin{1}: A figure handle to a figure to plot on to.
% varargin{2}: Plot type to use. Default is 'scatter'.

% Deal with optional inputs
if isempty(varargin)
    plotFig = figure;
else
    plotFig = varargin{1};
end

if (~isempty(varargin)) && (length(varargin) > 1)
    plotType = varargin{2};
else
    plotType = 'scatter';
end

% Specify plot variables
XVars(1).ProduceVar = @(st) st.SetSize;
XVars(1).NumBins = 'prebinned';

YVars(1).ProduceVar = @(st, inc) sum(st.Accuracy==1 & inc)/sum(inc);
YVars(1).FindIncludedTrials = @(st) true(size(st.DistractorMean));

YVars(2).ProduceVar = @(st, inc) sum(st.Response==1 & inc)/sum(inc);
YVars(2).FindIncludedTrials = @(st) st.Target == 1;

YVars(3).ProduceVar = @(st, inc) sum(st.Response==1 & inc)/sum(inc);
YVars(3).FindIncludedTrials = @(st) st.Target == 0;

Series(1).FindIncludedTrials = @(st) (st.BlockType==1);
Series(2).FindIncludedTrials = @(st) (st.BlockType==2);


PlotStyle.General = 'paper';

PlotStyle.Xaxis(1).Title = {'Number of items'};
PlotStyle.Xaxis(1).Ticks = [2, 3, 4, 6];

PlotStyle.Yaxis(1).Title = 'Accuracy';
PlotStyle.Yaxis(1).Ticks = linspace(0.4, 1, 7);
PlotStyle.Yaxis(1).TickLabels = string(PlotStyle.Yaxis(1).Ticks);
PlotStyle.Yaxis(1).InvisibleTickLablels = [2 : 2 : 7];
PlotStyle.Yaxis(1).RefVal = 0.5;

PlotStyle.Yaxis(2).Title = {'Hit rate'};
PlotStyle.Yaxis(2).Ticks = linspace(0.4, 1, 7);
PlotStyle.Yaxis(2).TickLabels = string(PlotStyle.Yaxis(2).Ticks);
PlotStyle.Yaxis(2).InvisibleTickLablels = [2 : 2 : 7];
PlotStyle.Yaxis(2).RefVal = 0.5;

PlotStyle.Yaxis(3).Title = {'FA rate'};
PlotStyle.Yaxis(3).Ticks = linspace(0, 1, 5);
PlotStyle.Yaxis(3).TickLabels = string(PlotStyle.Yaxis(3).Ticks);
PlotStyle.Yaxis(3).InvisibleTickLablels = [2 : 2 : 5];
PlotStyle.Yaxis(3).RefVal = 0.5;

PlotStyle.Data(1).Name = '{Uniform distractors}';
PlotStyle.Data(2).Name = '{Concentrated distractors}';
assert(all(unique(DSet.P(1).Data.KappaS( ...
    DSet.P(1).Data.BlockType==1))==0))

PlotStyle.Data(1).Colour = mT_pickColour(4);
PlotStyle.Data(2).Colour = mT_pickColour(2);


for i = 1 : length(PlotStyle.Data)
    PlotStyle.Data(i).PlotType = plotType;
end

plotFig = mT_plotVariableRelations(DSet, XVars, YVars, Series, ...
    PlotStyle, plotFig);
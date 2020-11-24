function plotFig = plotHitAndFaAllStats(DSet, varargin)

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


numBins = 10;

% Specify plot variables
XVars(1).ProduceVar = @(st) st.DistractorMean;
XVars(1).NumBins = numBins;
XVars(1).FindIncludedTrials = @(st) ~(st.SetSize==2); % Exclude trials
% for which there is only 2 items, because the
% variance of 1 distractor doesn't make sense

XVars(2).ProduceVar = @(st) st.DistractorVar;
XVars(2).NumBins = numBins;
XVars(2).FindIncludedTrials = @(st) ~(st.SetSize==2); % Exclude trials
% for which there is only 2 items, because the
% variance of 1 distractor doesn't make sense

XVars(3).ProduceVar = @(st) st.MostSimilarDistr;
XVars(3).NumBins = numBins;
XVars(3).FindIncludedTrials = @(st) ~(st.SetSize==2); % Exclude trials
% for which there is only 2 items, because the
% variance of 1 distractor doesn't make sense

YVars(1).ProduceVar = @(st, inc) sum(st.Response==1 & inc)/sum(inc);
YVars(1).FindIncludedTrials = @(st)  true(size(st.Target));

Series(1).FindIncludedTrials = @(st) (st.Target==1);
Series(2).FindIncludedTrials = @(st) (st.Target==0);


% Populate PlotStyle
PlotStyle.General = 'paper';

PlotStyle.Xaxis(1).Title = {'T-D mean', '(deg)'};
PlotStyle.Xaxis(1).Ticks = ...
    [0, pi/8, pi/4, 3*pi/8, pi/2, 5*pi/8, 3*pi/4, 7*pi/8, pi];
PlotStyle.Xaxis(1).TickLabels = {'0', ' ', ' ', ' ', '45', ' ', ' ', ' ', '90'};
PlotStyle.Xaxis(1).Lims = [-pi/16, pi];

PlotStyle.Xaxis(2).Title = 'D variance';
PlotStyle.Xaxis(2).Ticks = linspace(0, 1, 11);
PlotStyle.Xaxis(2).TickLabels = string(PlotStyle.Xaxis(2).Ticks);
PlotStyle.Xaxis(2).InvisibleTickLablels = setdiff([1:11], [1, 6, 11]);
PlotStyle.Xaxis(2).Lims = [-0.05, 1];

PlotStyle.Xaxis(3).Title = {'Min T-D', 'difference (deg)'};
PlotStyle.Xaxis(3).Ticks = ...
    [0, pi/8, pi/4, 3*pi/8, pi/2, 5*pi/8, 3*pi/4, 7*pi/8, pi];
PlotStyle.Xaxis(3).TickLabels = {'0', ' ', ' ', ' ', '45', ' ', ' ', ' ', '90'};
PlotStyle.Xaxis(3).Lims = [-pi/16, pi];


PlotStyle.Yaxis(1).Title = {'Probability ', '''present'' report'};
PlotStyle.Yaxis(1).Ticks = linspace(0, 1, 11);
PlotStyle.Yaxis(1).TickLabels = string(PlotStyle.Yaxis(1).Ticks);
PlotStyle.Yaxis(1).InvisibleTickLablels = setdiff(1:11, [1, 6, 11]);
PlotStyle.Yaxis(1).RefVal = 0.5;

PlotStyle.Data(1).Colour = [223, 172, 110]/255;
PlotStyle.Data(2).Colour = [132, 158, 209]/255;

PlotStyle.Data(1).PlotType = plotType;
PlotStyle.Data(2).PlotType = plotType;

PlotStyle.Data(1).Name = 'Target present (hit rate)';
PlotStyle.Data(2).Name = 'Target absent (FA rate)';

  
plotFig = mT_plotVariableRelations(DSet, XVars, YVars, Series, PlotStyle, plotFig);


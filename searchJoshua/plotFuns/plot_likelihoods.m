function Figures = plot_likelihoods(dataDir, parsDir)

[~, Nptpnts] = getData(dataDir);
Config = load('Config.mat');
Nmodels = length(Config.ModelList);

nLogLs = nan(Nptpnts,Nmodels);
stds = nan(Nptpnts,Nmodels);
for itype = 1 : Nmodels 
    fname = [parsDir '/pars_' Config.ModelList{itype} '.mat'];
    f = load(fname);

    sorted = sort(f.nLogLs,2);
    stds(:,itype) = std(sorted(:,1:3),[],2,'omitnan');
    nLogLs(:,itype) = min(f.nLogLs,[],2,'omitnan');
end

% First plot
pdata = nLogLs(:,1)-nLogLs;

xLocs = 0.5:2.5;
Figures.Likelihoods = figure;
bPlot = bar(xLocs, mean(pdata),'FaceColor','k')
hold on
plot(xLocs +0.05*randn(size(pdata)),pdata,'.','Color',[0.4,0.4,0.4])
box off
set(gca,'TickDir','out')
ylabel('LL - LLBayes')
set(gca)
xlim([0, 3])
xticks(xLocs)

% Set labels, although replace some to make more interpretable
xticklabels(Config.ModelLabel);

set(findall(gcf,'-property','FontSize'),'FontSize',10)
set(findall(gcf,'Type','Line'), 'LineWidth', 1)
set(findall(gca,'Type','Line'), 'LineWidth', 1)
bPlot.BaseLine.LineWidth = 1;
set(gca, 'LineWidth', 1)


% Second plot
pdata2 = stds;

figure 
bar(mean(pdata2),'FaceColor','k')
hold on
plot((1:3) +0.05*randn(size(pdata2)),pdata2,'.','Color',[0.4,0.4,0.4],'MarkerSize',20)
box off
set(gca,'TickDir','out')
ylabel('std of estimates','FontSize',18)
set(gca,'FontSize',14)
xticklabels(Config.ModelList)


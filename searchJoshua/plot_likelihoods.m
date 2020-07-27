function Figures = plot_likelihoods(dataDir)

[~, Nptpnts] = getData(dataDir);
Config = load('Config.mat');
Nmodels = length(Config.ModelList);

nLogLs = nan(Nptpnts,Nmodels);
stds = nan(Nptpnts,Nmodels);
for itype = 1 : Nmodels 
    fname = ['pars/pars_' Config.ModelList{itype} '.mat'];
    f = load(fname);

    sorted = sort(f.nLogLs,2);
    stds(:,itype) = std(sorted(:,1:3),[],2,'omitnan');
    nLogLs(:,itype) = min(f.nLogLs,[],2,'omitnan');
end

pdata = nLogLs(:,1)-nLogLs;

Figures.Likelihoods = figure;
bar(mean(pdata),'FaceColor','k')
hold on
plot((1:3) +0.05*randn(size(pdata)),pdata,'.','Color',[0.4,0.4,0.4],'MarkerSize',20)
box off
set(gca,'TickDir','out')
ylabel('LL - LLBayes','FontSize',18)
set(gca,'FontSize',14)
xticklabels(Config.ModelList)


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
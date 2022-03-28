function plot_likelihoods(Naverage, subset)

if ~exist('subset','var') || isempty(subset)
    subset = 1:4;
end

Nsubs = 9;

likelihoods = zeros(Nsubs,4);
stds = zeros(Nsubs,4);
for itype = 1:4
    switch itype
        case 1
            fname = 'pars/pars_ibs_Bayes.mat';
        case 2
            fname = 'pars/pars_ibs_Freq.mat';
        case 3
            fname = 'pars/pars_ibs_Sample.mat';
        case 4
            fname = 'pars/pars_ibs_cssample.mat';
    end
    f = load(fname);
    if itype>5
        sorted = sort(f.likelihoods,2);
        if exist('Naverage','var') && ~isempty(Naverage)
            stds(:,itype) = std(sorted(:,1:Naverage),[],2,'omitnan')/log(2);
            likelihoods(:,itype) = mean(sorted(:,1:Naverage),2,'omitnan')/log(2);
        else
            stds(:,itype) = std(sorted(:,1:3),[],2,'omitnan')/log(2);
            likelihoods(:,itype) = min(f.likelihoods,[],2,'omitnan')/log(2);
        end
    else
        sorted = sort(f.likelihoods,2);
        if exist('Naverage','var') && ~isempty(Naverage)
            stds(:,itype) = std(sorted(:,1:Naverage),[],2,'omitnan');
            likelihoods(:,itype) = mean(sorted(:,1:Naverage),2,'omitnan');
        else
            stds(:,itype) = std(sorted(:,1:3),[],2,'omitnan');
            likelihoods(:,itype) = min(f.likelihoods,[],2,'omitnan');
        end
    end
end

pdata = likelihoods(:,1)-likelihoods;

figure 
bar(mean(pdata(:, subset)),'FaceColor','k')
hold on
plot((1:length(subset)) +0.05*randn(size(pdata(:,subset))),pdata(:,subset),'.','Color',[0.4,0.4,0.4],'MarkerSize',20)
box off
set(gca,'TickDir','out', 'LineWidth', 2)
ylabel('LL - LLBayes','FontSize',20)
set(gca,'FontSize', 18, 'XColor', 'k', 'YColor', 'k')
labels = {'Bayes', 'Point Estimate', 'Sampling', 'joint Sampling'};
xticklabels(labels(subset));
xlim([0.5, length(subset) + 0.5])
set(gcf, 'Position', [550   750   560   250])


pdata2 = stds;

figure 
bar(mean(pdata2),'FaceColor','k')
hold on
plot((1:length(subset)) +0.05*randn(size(pdata2(:,subset))),pdata2(:,subset),'.','Color',[0.4,0.4,0.4],'MarkerSize',20)
box off
set(gca,'TickDir','out')
ylabel('std of estimates','FontSize',20)
set(gca,'FontSize',14)
xticklabels({'Bayes','Point Estimate', 'Sampling', 'joint Sampling'})
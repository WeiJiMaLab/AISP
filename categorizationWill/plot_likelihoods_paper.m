function plot_likelihoods_paper(Naverage)

Nsubs = 11;

likelihoods = zeros(Nsubs,3);
stds = zeros(Nsubs,3);
for itype = 1:3
    switch itype
        case 1
            fname = 'pars/pars_ibs_Bayes.mat';
        case 2
            fname = 'pars/pars_ibs_Freq.mat';
        case 3
            fname = 'pars/pars_ibs_Var.mat';
    end
    f = load(fname);
    if itype>5
        sorted = sort(f.likelihoods,2);
        if exist('Naverage','var')
            stds(:,itype) = std(sorted(:,1:Naverage),[],2,'omitnan')/log(2);
            likelihoods(:,itype) = mean(sorted(:,1:Naverage),2,'omitnan')/log(2);
        else
            stds(:,itype) = std(sorted(:,1:3),[],2,'omitnan')/log(2);
            likelihoods(:,itype) = min(f.likelihoods,[],2,'omitnan')/log(2);
        end
    else
        sorted = sort(f.likelihoods,2);
        if exist('Naverage','var')
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
bar(mean(pdata),'FaceColor','k')
hold on
plot((1:3) +0.05*randn(size(pdata)),pdata,'.','Color',[0.4,0.4,0.4],'MarkerSize',20)
box off
set(gca,'TickDir','out')
ylabel('LL - LLBayes','FontSize',18)
set(gca,'FontSize', 14)
xticklabels({'Bayes','PE','Vari'})


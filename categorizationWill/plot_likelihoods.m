function plot_likelihoods

Nsubs = 11;

likelihoods = zeros(Nsubs,5);
for itype = 1:5
    switch itype
        case 1
            fname = 'pars/parsBayes.mat';
        case 2
            fname = 'pars/parsFreq.mat';
        case 3
            fname = 'pars/parsFreq2.mat';
        case 4
            fname = 'pars/parsFreq3.mat';
        case 5
            fname = 'pars/parsVar.mat';
    end
    f = load(fname);
    likelihoods(:,itype) = min(f.likelihoods,[],2);
end

pdata = likelihoods(:,1)-likelihoods;

figure 
bar(mean(pdata),'FaceColor','k')
hold on
plot((1:5) +0.05*randn(size(pdata)),pdata,'.','Color',[0.4,0.4,0.4],'MarkerSize',20)
box off
set(gca,'TickDir','out')
ylabel('LL - LLBayes','FontSize',18)
set(gca,'FontSize',14)
xticklabels({'Bayes','PE1','PE2','PE3','Vari'})
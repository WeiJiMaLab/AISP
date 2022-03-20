function plot_check(slack, N_target)

if ~exist('slack', 'var') || isempty(slack)
    slack = 2;
end
if ~exist('N_target', 'var') || isempty(N_target)
    N_target = 5;
end

Nsubs = 9;

likelihoods = zeros(Nsubs,5);
highest = zeros(Nsubs,5);
n_ok = zeros(Nsubs,5);
for itype = 1:5
    switch itype
        case 1
            fname = 'pars/pars_ibs_Bayes.mat';
        case 2
            fname = 'pars/pars_ibs_Freq.mat';
        case 4
            fname = 'pars/pars_ibs_Sample.mat';
        case 5
            fname = 'pars/pars_ibs_cssample.mat';
    end
    f = load(fname);
    sorted = sort(f.likelihoods,2);
    n_ok(:,itype) = nansum((sorted-sorted(:,1)) < slack, 2);
    highest(:,itype) = sorted(:, N_target)- sorted(:, 1);
    likelihoods(:,itype) = min(f.likelihoods,[],2,'omitnan');
end

figure
imagesc(highest)
colorbar()
box off
set(gca,'TickDir','out')
ylabel('Subject','FontSize',18)
set(gca,'FontSize', 14)
xticks(1:5)
xticklabels({'Bayes','PE','PE2', 'Sample', 'cssample'})

figure
imagesc(n_ok,[1,5])
colorbar()
box off
set(gca,'TickDir','out')
ylabel('Subject','FontSize',18)
set(gca,'FontSize', 14)
xticks(1:5)
xticklabels({'Bayes','PE','PE2', 'Sample', 'cssample'})

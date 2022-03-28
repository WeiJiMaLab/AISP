function plot_performance(N, beta_1)
% plots the performance in % correct against the amount of observer noise

if ~exist('N','var') || isempty(N)
    N = 1000;
end

pars_B = load('pars/pars_ibs_Bayes.mat');
pars = pars_B.pars;
lik = pars_B.likelihoods;
[l_min, idx] = nanmin(lik, [], 2);
p_best = pars(:, :, 1);
for i = 1:size(pars, 1)
    p_best(i, :) = pars(i, :, idx(i));
end

sigmas = exp(mean(p_best(:,1:6),1));
beta = mean(pars(:,7:8));
if ~exist('beta_1','var') || isempty(beta_1)
    beta(2) = exp(beta(2));
else
    beta(2) = beta_1;
end
beta(1) = 0;
lambda = 0;
sigmas_test = exp(linspace(log(min(sigmas)), max(log(sigmas)), 100));

sigma1 = 3;
sigma2 = 12;
data = [ ...
    [ones(N,1), sigma1 .* randn(N,1)];
    [ones(N,1), sigma2 .* randn(N,1)];
    ];
sigmas_t = sigmas;
performances = zeros(length(sigmas_test), 5);
cat_true = [zeros(N,1); ones(N,1)];
k = 1;
for i_sigma = sigmas_test
    sigmas_t(1) = i_sigma;
    response_bayes = bayes_simulate(data, sigmas_t, beta, lambda);
    response_freq = freq_simulate(data, sigmas_t, beta, lambda);
    response_var = var_simulate(data, sigmas_t, beta, lambda);
    response_samp = sample_simulate(data, sigmas_t, beta, lambda, 10);
    response_cssamp = cssample_simulate(data, sigmas_t, beta, lambda, 10);
    performances(k, :) = [ ...
        mean(cat_true == response_bayes),...
        mean(cat_true == response_freq),...
        mean(cat_true == response_var),...
        mean(cat_true == response_samp),...
        mean(cat_true == response_cssamp)];
    k = k + 1;
end

figure
hold on
colors = [0  ,0  ,0;...
     0  ,0  ,1;...
     0  ,1  ,0;...
     1  ,0  ,0;...
     1  ,0.7,0.7];
for i = 1:5
    plot(sigmas_test, performances(:,i), 'LineWidth', 3, 'Color', colors(i,:))
end
set(gca,'xscale' , 'log')
yline(0.5, '--', 'LineWidth', 2)
for s = sigmas
    plot([s,s], [0.45,0.5], 'k-', 'LineWidth', 2)
end
ylim([0.4,1])
xlim([0.3,40])
xlabel('\sigma noise', 'FontSize', 18)
ylabel('Proportion correct', 'FontSize', 18)
legend({'Bayes', 'point estimate', 'variational', 'sample 10', 'joint sample 10'})
legend boxoff
box off
set(gca, 'TickDir', 'out', 'FontSize', 16, 'LineWidth', 2, 'GridColor', 'k');
set(gca, 'box', 'off');
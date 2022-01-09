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
    response_bayes = bayes_simulate_vec(data, sigmas_t, beta, lambda);
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
semilogx(sigmas_test, performances, 'LineWidth', 2)
hold on
plot([0.1, 100], [0.5,0.5], 'k--', 'LineWidth', 2)
for s = sigmas
    plot([s,s], [0.475,0.5], 'k-', 'LineWidth', 1)
end
set(gca, 'FontSize', 16, 'TickDir', 'out')
xlabel('\sigma noise', 'FontSize', 20)
ylabel('Proportion Correct', 'FontSize', 20)
legend({'Bayes', 'point estimate', 'variational', 'sample 10', 'joint sample 10'})
legend boxoff
box off
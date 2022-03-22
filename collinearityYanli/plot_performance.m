function plot_performance(N, beta_1, n_sample)
% plots the performance in % correct against the amount of observer noise

if ~exist('N','var') || isempty(N)
    N = 1000;
end
if ~exist('n_sample','var') || isempty(n_sample)
    n_sample = 10;
end

pars_B = load('pars/pars_Bayes.mat');
pars = pars_B.pars;
lik = pars_B.likelihoods;
[l_min, idx] = nanmin(lik, [], 2);
p_best = pars(:, :, 1);
for i = 1:size(pars, 1)
    p_best(i, :) = pars(i, :, idx(i));
end

sigmas = mean(p_best(:, 1:4), 1);
beta = mean(p_best(:,5:6), 1);
if ~exist('beta_1','var') || isempty(beta_1)
    beta(2) = exp(beta(2));
else
    beta(2) = beta_1;
end
beta(1) = 0;
lambda = 0;
sigmas_test = exp(linspace(log(0.25), log(50), 100));

sigma0 = 24;
data = [ ...
    [zeros(N,1), sigma0 .* randn(N,2)];
    [zeros(N,1), repmat(sigma0 .* randn(N,1),1,2)];
    ];
sigmas_t = sigmas;
performances = zeros(length(sigmas_test), 5);
cat_true = [zeros(N,1); ones(N,1)];
k = 1;
for i_sigma = sigmas_test
    sigmas_t(1) = i_sigma;
    response_bayes = bayes_simulate(data, [sigmas_t, beta, lambda]);
    response_freq = pe_simulate(data, [sigmas_t, beta, lambda]);
    response_freq2 = pe2_simulate(data, [sigmas_t, beta, lambda]);
    response_samp = sample_simulate(data, [sigmas_t, beta, lambda, n_sample]);
    response_cssamp = cssample_simulate(data, [sigmas_t, beta, lambda, n_sample]);
    performances(k, :) = [ ...
        mean(cat_true == response_bayes),...
        mean(cat_true == response_freq),...
        mean(cat_true == response_freq2),...
        mean(cat_true == response_samp),...
        mean(cat_true == response_cssamp)];
    k = k + 1;
end

figure
hold on
colors = [0  ,0  ,0;...
     0.7,0.7,1;...
     0  ,0  ,1;...
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
set(gca, 'FontSize', 16, 'TickDir', 'out', 'LineWidth', 2, 'GridColor', 'k')
xlabel('\sigma noise', 'FontSize', 18)
ylabel('Proportion Correct', 'FontSize', 18)
legend({'Bayes', 'point estimate', 'point estimate opt.', ...
    sprintf('sample %d', n_sample), sprintf('joint sample %d', n_sample)})
legend boxoff
box off
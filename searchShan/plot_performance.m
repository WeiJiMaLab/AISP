function plot_performance(N, beta_1, n_sample)
% plots the performance in % correct against the amount of observer noise

if ~exist('N','var') || isempty(N)
    N = 1000;
end

if ~exist('n_sample','var') || isempty(n_sample)
    n_sample = 10;
end

pars_B = load('pars/pars_ibs_Bayes.mat');
pars = pars_B.pars;
lik = pars_B.likelihoods;
[l_min, idx] = nanmin(lik, [], 2);
p_best = pars(:, :, 1);
for i = 1:size(pars, 1)
    p_best(i, :) = pars(i, :, idx(i));
end

sigmas = mean(p_best(:, 1), 1);
beta = mean(p_best(:,2:3), 1);
if ~exist('beta_1','var') || isempty(beta_1)
    beta(1) = beta(1);
else
    beta(1) = beta_1;
end
beta(2) = 0;
lambda = 0;
sigmas_test = exp(linspace(log(0.25), log(50), 100));

sigma = 9.06;
data = sigma .* randn(2*N,2);
data(1:N,1) = -abs(data(1:N,1));
data((N+1):(2*N),1) = abs(data((N+1):(2*N),1));
performances = zeros(length(sigmas_test), 4);
cat_true = [zeros(N,1); ones(N,1)];
k = 1;
for i_sigma = sigmas_test
    response_bayes = bayes_simulate(data, [i_sigma, beta, lambda]);
    response_freq = pe_simulate(data, [i_sigma, beta, lambda]);
    response_samp = sample_simulate(data, [i_sigma, beta, lambda, n_sample]);
    response_cssamp = cssample_simulate(data, [i_sigma, beta, lambda, n_sample]);
    performances(k, :) = [ ...
        mean(cat_true == response_bayes),...
        mean(cat_true == response_freq),...
        mean(cat_true == response_samp),...
        mean(cat_true == response_cssamp)];
    k = k + 1;
end

figure
hold on
colors = [0  ,0  ,0;...
     0  ,0  ,1;...
     1  ,0  ,0;...
     1  ,0.7,0.7];
for i = 1:4
    plot(sigmas_test, performances(:,i), 'LineWidth', 3, 'Color', colors(i,:))
end
set(gca,'xscale' , 'log')
yline(0.5, '--', 'LineWidth', 2)
for s = sigmas
    plot([s,s], [0.45,0.5], 'k-', 'LineWidth', 2)
end
set(gca, 'FontSize', 16, 'TickDir', 'out', 'LineWidth', 2, 'GridColor', 'k')
xlabel('\sigma noise', 'FontSize', 18)
ylabel('Proportion Correct', 'FontSize', 18)
legend({'Bayes', 'point estimate', ...
    sprintf('sample %d', n_sample), sprintf('joint sample %d', n_sample)})
legend boxoff
box off
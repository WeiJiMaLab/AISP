function plot_performance(N, beta_1, n_sample)
% plots the performance in % correct against the amount of observer noise

addpath(genpath('helper_functions'))

if ~exist('N','var') || isempty(N)
    N = 10000;
end

if ~exist('n_sample','var') || isempty(n_sample)
    n_sample = 10;
end

tempp = load('cdf_table.mat');

pars_B = load('fits/modelbayes.mat');
pars = pars_B.bfpMat;
pars = mean(pars);
sigma_fit = pars(2);

if ~exist('beta_1','var') || isempty(beta_1)
    pars(4) = pars(4);
else
    pars(4) = beta_1;
end
% bias and lapse = 0
pars(5) = 0;
pars(6) = 0;
pars(7) = n_sample;

prec_factor = pars(1) / pars(2);
prec_test = exp(linspace(log(1), log(200), 20));

nItems = 4;
DeltaMat = zeros(nItems, 2*N);
DeltaVec = (rand(N, 1) .* 2 .* pi) - pi; % all Deltas
idx = randi(4, [N, 1]) ...
    + (0:4:(N * nItems - nItems))';
DeltaMat(idx) = DeltaVec(:);
DeltaMat = DeltaMat';

rel = randi(2, 2*N, 4);

data = [DeltaMat, rel];

cat_true = [ones(N,1); zeros(N,1)];
k = 1;
performances = zeros(length(prec_test), 5);
for i_prec = prec_test
    pars(1) = prec_factor * i_prec;
    pars(2) = i_prec;
    response_bayes = simulate_responses(pars,'bayes',data,[],tempp);
    % Decision variable for frequentist is always > 0, some bias required
    response_freq = simulate_responses(pars + [0,0,0,0,10,0,0],'freq',data,[],tempp);
    response_freq2 = simulate_responses(pars,'freq2',data,[],tempp);
    response_samp = sample_simulate2(pars,data,[],tempp);
    response_cssamp = cssample_simulate3(pars,data,[],tempp);
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
    plot(prec_test, performances(:,i), 'LineWidth', 3, 'Color', colors(i,:))
end
set(gca,'xscale' , 'log', 'XDir', 'reverse')
yline(0.5, '--', 'LineWidth', 2)
plot([sigma_fit,sigma_fit], [0.45,0.5], 'k-', 'LineWidth', 3);
set(gca, 'FontSize', 16, 'TickDir', 'out', 'LineWidth', 2, 'GridColor', 'k')
xlabel('Average precision \kappa', 'FontSize', 18)
ylabel('Proportion Correct', 'FontSize', 18)
xlim([0.7,600])
legend({'Bayes', 'point estimate', 'point estimate opt.' ...
    sprintf('sample %d', n_sample), sprintf('joint sample %d', n_sample)})
legend boxoff
box off
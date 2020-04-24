function plot_decisions(type, pars, N)

if ~exist('N', 'var')|| isempty(N)
    N = 100000;
end

if ~exist('pars', 'var') || isempty(pars)
    pars = [5,20,0,0]; % idealized 
end

N_quantiles = 9;

sigma = 9.06;

stimulus = sigma * randn(N,2);

switch type
    case 'Bayes'
        response = bayes_simulate(stimulus, pars);
    case 'Freq'
        response = pe_simulate(stimulus, pars);
end

bin = 1/(N_quantiles);
bins = [-inf,9.06 * norminv(linspace(bin,1-bin,N_quantiles-1)),inf];


count_true = zeros(N_quantiles, N_quantiles);
count_all = zeros(N_quantiles, N_quantiles);
for i_bin = 1:N_quantiles
    for j_bin = 1:N_quantiles
        select = (stimulus(:, 1) > bins(i_bin)) & (stimulus(:, 1) <= bins(i_bin + 1)) ...
            & (stimulus(:, 2) > bins(j_bin)) & (stimulus(:, 2) <= bins(j_bin + 1));
        count_true(i_bin, j_bin) = sum(response(select));
        count_all(i_bin, j_bin) = sum(select);
    end
end

imagesc((count_true./count_all)', [0,1])
axis square
set(gca,'XTick',[])
set(gca,'YTick',[])
colorbar()

function plot_decisions_lines(type)

switch type
    case 'Bayes'
        pars = load('pars/pars_ibs_Bayes.mat');
        lik = pars.likelihoods;
        pars = pars.pars;
        [~, idx] = min(lik,[],2);
        pars_sub = nan(size(pars,1),size(pars,2));
        for i_sub = 1:size(pars,1)
        	pars_sub(i_sub,:) = pars(i_sub,:,idx(i_sub));
        end
    case 'Freq'
        pars = load('pars/pars_ibs_Freq.mat');
        lik = pars.likelihoods;
        pars = pars.pars;
        [~, idx] = min(lik,[],2);
        pars_sub = nan(size(pars,1),size(pars,2));
        for i_sub = 1:size(pars,1)
        	pars_sub(i_sub,:) = pars(i_sub,:,idx(i_sub));
        end
end

% Colors from aspen
purple = [150 112 177] / 255; %#9670B1
gold = [223 172 110] / 255; %#DFAC6E
reddish = [213 96 100] / 255; %#D56064
blue = [132 158 209] / 255; %#849ED1

N_quantiles_y = 9;
N_quantiles_x = 5;

colors = gold + linspace(0, 1, N_quantiles_x)' * (blue - gold);


sigma = 9.06;

response_data = [];
stimulus_data = [];
response = nan(9,3000);
for i_sub = 1:9
    [stim, resp, perf] = readdata(i_sub);
    response_data = cat(2, response_data, resp);
    stimulus_data = cat(3, stimulus_data, stim);
    switch type
        case 'Bayes'
            response(i_sub, :) = bayes_simulate(stim, pars_sub(i_sub,:));
        case 'Freq'
            response(i_sub, :) = pe_simulate(stim, pars_sub(i_sub,:));
    end
end
response_data = 0.5 + 0.5 * response_data;


bin = 1/(N_quantiles_x);
bins_x = [-inf,9.06 * norminv(linspace(bin,1-bin,N_quantiles_x-1)),inf];
bin = 1/(N_quantiles_y);
bins_y = [-inf,9.06 * norminv(linspace(bin,1-bin,N_quantiles_y-1)),inf];


count_true = zeros(N_quantiles_x, N_quantiles_y, 9);
count_all = zeros(N_quantiles_x, N_quantiles_y, 9);
count_true_data = zeros(N_quantiles_x, N_quantiles_y, 9);
count_all_data = zeros(N_quantiles_x, N_quantiles_y, 9);
for i_bin = 1:N_quantiles_x
    for j_bin = 1:N_quantiles_y
        select_data = (stimulus_data(:, 1, :) > bins_x(i_bin)) & (stimulus_data(:, 1, :) <= bins_x(i_bin + 1)) ...
            & (stimulus_data(:, 2, :) > bins_y(j_bin)) & (stimulus_data(:, 2, :) <= bins_y(j_bin + 1));
        for i_sub = 1:9
            response_sub = response(i_sub,:);
            count_true(i_bin, j_bin, i_sub) = sum(response_sub(select_data(:,1,i_sub)));
            count_all(i_bin, j_bin, i_sub) = sum(select_data(:,:,i_sub));
            resp_sub = response_data(:, i_sub);
            count_true_data(i_bin, j_bin, i_sub) = sum(resp_sub(select_data(:,1,i_sub)));
            count_all_data(i_bin, j_bin, i_sub) = sum(select_data(:,:,i_sub));
        end
    end
end

prop_data = count_true_data ./ count_all_data;
prop_mean_data = mean(prop_data, 3,'omitnan');
prop_std_data = std(prop_data, [], 3,'omitnan')./sqrt(9);
prop = count_true ./ count_all;
prop_mean = mean(prop, 3,'omitnan');
prop_std = std(prop, [], 3,'omitnan')./sqrt(9);

figure
hold on
for i_quantile = 1:N_quantiles_x
    col = colors(i_quantile, :);
    s_bin = 1:N_quantiles_y;
    x = [s_bin';flip(s_bin)'];
    y = [prop_mean(i_quantile,:)-prop_std(i_quantile,:), ...
        flip(prop_mean(i_quantile,:)+prop_std(i_quantile,:))]';
    plot(polyshape(x,y), 'LineStyle', 'none', 'FaceColor', col, 'FaceAlpha', 0.5)
    errorbar(1:N_quantiles_y, prop_mean_data(i_quantile,:), prop_std_data(i_quantile,:), '.', 'Color', col, ...
      'MarkerSize', 10, 'LineWidth', 2)
end

hleg = legend({'far left','','left','','near vertical','','right','','far right',''}, 'box','off', 'FontSize', 16);
title(hleg,'Target orientation bins', 'FontSize', 14);
ylim([0,1])
yax = get(gca, 'YAxis');
set(yax, 'TickLength', [0.02,0.05])
set(yax, 'TickValues', [0,0.5,1])
set(yax, 'MinorTick', 'on')
set(yax, 'MinorTickValues', [0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9])
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 2)
xax = get(gca, 'XAxis');
set(xax, 'TickLength', [0.02,0.05])
set(xax, 'TickValues', [1,5,9])
set(xax, 'MinorTick', 'on')
set(xax, 'MinorTickValues', 1:9)
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 2)
xlabel('Distractor orientation [bins]', 'FontSize',20)
ylabel('Proportion reports "right"', 'FontSize',20)
set(hleg, 'Position', [0.75,0.25,0.25,0.5])
set(gcf, 'Position', [560, 530, 650, 420])
set(gca, 'Position', [0.1300 0.14 0.6 0.8])

function plot_offset(data, n_bin, model, n_repeat)

if ~exist('n_bin', 'var') || isempty(n_bin)
    n_bin = 9;
end

if ~exist('model', 'var') || isempty(model)
    model = 'None';
end

if ~exist('n_repeat', 'var') || isempty(n_repeat)
    n_repeat = 25;
end
% Colors from aspen
purple = [150 112 177] / 255; %#9670B1
gold = [223 172 110] / 255; %#DFAC6E
reddish = [213 96 100] / 255; %#D56064
blue = [132 158 209] / 255; %#849ED1
colors = blue + linspace(0, 1, 4)' * (gold - blue);

bins = [-inf, linspace(-50, 50, n_bin-1), inf];
subs = unique(data(:,1));

figure
hold on

data_model = data;
if ~strcmp(model, 'None')
    if strcmp(model, 'bayes')
        p_data = load('pars/pars_Bayes.mat');
    elseif strcmp(model, 'PE')
        p_data = load('pars/pars_Freq.mat');
    elseif strcmp(model, 'PE2')
        p_data = load('pars/pars_Freq2.mat');
    end
    pars = p_data.pars;
    lik = p_data.likelihoods;
    p_sub_model = zeros(length(subs), n_bin, 4);
    for i_sub = subs'
        [val, idx] = min(lik(i_sub,:));
        responses = zeros(sum(data_model(:,1)==i_sub), n_repeat);
        dat_sub = data_model(data_model(:,1)==i_sub, :);
        for i_rep = 1:n_repeat
            if strcmp(model, 'bayes')
                responses(:,i_rep) = bayes_simulate( ...
                    dat_sub(:,[2, 5, 6]), pars(i_sub,:,idx));
            elseif strcmp(model, 'PE')
                responses(:,i_rep) = pe_simulate( ...
                    dat_sub(:,[2, 5, 6]), pars(i_sub,:,idx));
            elseif strcmp(model, 'PE2')
                responses(:,i_rep) = pe2_simulate( ...
                    dat_sub(:,[2, 5, 6]), pars(i_sub,:,idx));
            end
        end
        data_model(data_model(:,1)==i_sub, 4) = mean(responses, 2);
        sdiff = dat_sub(:,6) - dat_sub(:,5);
        for i_bin = 1:n_bin
            idx_bin = (sdiff > bins(i_bin)) & (sdiff < bins(i_bin+1));
            p_sub_model(i_sub, i_bin, 1) = mean(mean(responses(idx_bin & (dat_sub(:,2)==0),:)));
            p_sub_model(i_sub, i_bin, 2) = mean(mean(responses(idx_bin & (dat_sub(:,2)==240),:)));
            p_sub_model(i_sub, i_bin, 3) = mean(mean(responses(idx_bin & (dat_sub(:,2)==480),:)));
            p_sub_model(i_sub, i_bin, 4) = mean(mean(responses(idx_bin & (dat_sub(:,2)==840),:)));
        end
    end
end


for i_ecc = 1:4
    col = colors(i_ecc, :);
    switch i_ecc
        case 1
            dat = data(data(:,2)==0,:);
        case 2
            dat = data(data(:,2)==240,:);
        case 3
            dat = data(data(:,2)==480,:);
        case 4
            dat = data(data(:,2)==840,:);
    end
    %s = mean(dat(:,5:6),2) - dat(:,2);
    sdiff = dat(:,6) - dat(:,5);
    s_bin = zeros(n_bin, 1);
    p_same = zeros(n_bin, 1);
    p_std = zeros(n_bin, 1);
    for i_bin = 1:n_bin
        dat_bin = dat((sdiff > bins(i_bin)) & (sdiff < bins(i_bin+1)), :);
        p_same(i_bin) = mean(dat_bin(:,4));
        s_bin(i_bin) = mean(sdiff((sdiff > bins(i_bin)) & (sdiff < bins(i_bin+1))));
        % calculate the standard deviation
        p_sub = zeros(length(subs), 1);
        k_sub = 0;
        for i_sub = subs'
            k_sub = k_sub + 1;
            p_sub(k_sub) = mean(dat_bin(dat_bin(:,1)==i_sub,4));
        end
        p_std(i_bin) = std(p_sub);
    end
    if ~strcmp(model, 'None')
        p_model = mean(p_sub_model(:,:,i_ecc),1);
        std_model = std(p_sub_model(:,:,i_ecc),1);
        x = [s_bin;flip(s_bin)];
        y = [p_model-std_model./sqrt(length(subs)),flip(p_model+std_model./sqrt(length(subs)))]';
        plot(polyshape(x,y), 'LineStyle', 'none', 'FaceColor', col, 'FaceAlpha', 0.5)
    end
    errorbar(s_bin, p_same, p_std ./ sqrt(length(subs)), '.', 'Color', col, ...
        'MarkerSize', 10, 'LineWidth', 2)
end

hleg = legend({'0','','240','','480','','840',''}, 'box','off', 'FontSize', 16);
title(hleg,'Eccentricity[pix]', 'FontSize', 14);
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
set(xax, 'TickValues', [-100,-50,0,50,100])
set(xax, 'MinorTick', 'on')
set(xax, 'MinorTickValues', -100:10:100)
set(gca, 'TickDir', 'out')
set(gca, 'FontSize', 16)
set(gca, 'LineWidth', 2)
xlabel('Offset [pixel]', 'FontSize',20)
ylabel('Proportion same reports', 'FontSize',20)
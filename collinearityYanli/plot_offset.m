function plot_offset(data, n_bin)

if ~exist('n_bin', 'var') || isempty(n_bin)
    n_bin = 9;
end

bins = [-inf, linspace(-50, 50, n_bin-1), inf];

figure
hold on
for i_ecc = 1:4
    switch i_ecc
        case 1
            dat = data(data(:,2)==0,:);
            col = 'r';
        case 2
            dat = data(data(:,2)==240,:);
            col = 'g';
        case 3
            dat = data(data(:,2)==480,:);
            col = 'b';
        case 4
            dat = data(data(:,2)==840,:);
            col = 'k';
    end
    %s = mean(dat(:,5:6),2) - dat(:,2);
    sdiff = dat(:,6) - dat(:,5);
    s_bin = zeros(n_bin, 1);
    p_same = zeros(n_bin, 1);
    for i_bin = 1:n_bin
        dat_bin = dat((sdiff > bins(i_bin)) & (sdiff < bins(i_bin+1)), :);
        p_same(i_bin) = mean(dat_bin(:,4));
        s_bin(i_bin) = mean(sdiff((sdiff > bins(i_bin)) & (sdiff < bins(i_bin+1))));
    end
    plot(s_bin, p_same, '.-', 'Color', col, 'MarkerSize', 10, 'LineWidth', 2);
end

set(gca, 'TickDir', 'out')
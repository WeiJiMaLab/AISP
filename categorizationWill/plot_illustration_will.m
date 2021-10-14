function plot_illustration_will
clf

sigma1 = 3;
sigma2 = 12;

x = linspace(-30,30,1000);
y1 = exp(-0.5 * (x).^2 ./ sigma1.^2) ./ sqrt(2 *pi) ./ sigma1;
y2 = exp(-0.5 * (x).^2 ./ sigma2.^2) ./ sqrt(2 *pi) ./ sigma2;

plot(x, y1, 'LineWidth', 2);
hold on
plot(x, y2, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'Fontsize', 14, 'LineWidth',2)
set(gca, 'XTick', [-30, -20, -10, 0, 10, 20, 30])
set(gca, 'YTick', [])
xlabel('Stimulus Level s')
ylabel('Density')
legend('$C=0$', '$C=1$', 'Location', 'northwest', ...
    'Interpreter','latex')
legend('boxoff')
box off
xlim([-30,30])
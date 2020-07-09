function plot_illustration

sigma1 = 3;
sigma2 = 12;

x = linspace(-30,30,1000);
y1 = exp(-0.5 * x.^2 ./ sigma1.^2) ./ sqrt(2 *pi) ./ sigma1;
y2 = exp(-0.5 * x.^2 ./ sigma2.^2) ./ sqrt(2 *pi) ./ sigma2;

hold on
plot(x, y1, 'LineWidth', 2);
plot(x, y2, 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'Fontsize', 14)
xlabel('Stimulus Level')
function plot_illustration
clf

sigma1 = 6;
sigma2 = 12;
sigmaL = 3;
x0 = 3;

x = linspace(-30,30,1000);
y1 = exp(-0.5 * (x+10).^2 ./ sigma1.^2) ./ sqrt(2 *pi) ./ sigma1;
y2 = exp(-0.5 * (x-10).^2 ./ sigma1.^2) ./ sqrt(2 *pi) ./ sigma1;
l = exp(-0.5 * (x-x0).^2 ./ sigmaL.^2) ./ sqrt(2 *pi) ./ sigmaL;

subplot(2,3,1)
hold on
plot(x, y1, 'LineWidth', 2);
plot(x, y2, 'LineWidth', 2);
plot(x, l, '--k', 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'Fontsize', 14)
set(gca, 'XTick', [0, x0])
set(gca, 'YTick', [])
set(gca, 'XTickLabel', {'0', '$x$'}, 'TickLabelInterpreter','latex')
xlabel('Stimulus Level s')
ylabel('Density')
legend('$C=0$', '$C=1$', 'likelihood', 'Location', 'northwest', ...
    'Interpreter','latex')
legend('boxoff')
xlim([-25,25])
title('Prior & likelihood')

subplot(2,3,2)
hold on
patch(x, y2 .* l, 'red', 'FaceAlpha', 0.3, 'LineWidth', 2, 'EdgeColor', 'red');
patch(x, y1 .* l, 'blue', 'FaceAlpha', 0.3, 'LineWidth', 2, 'EdgeColor', 'blue');
[m1, x1] = max(y1 .* l);
[m2, x2] = max(y2 .* l);
x1 = x(x1);
x2 = x(x2);
plot([-30,x1], [m1, m1], 'k:', 'LineWidth', 2);
plot([-30,x2], [m2, m2], 'k:', 'LineWidth', 2);
plot([x1,x1], [0, m1], 'k:', 'LineWidth', 2);
plot([x2,x2], [0, m2], 'k:', 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'Fontsize', 14)
set(gca, 'XTick', [x1, x0, x2])
set(gca, 'YTick', [])
set(gca, 'XTickLabel', {'$\hat{s}_0$','', '$\hat{s}_1$'}, 'TickLabelInterpreter','latex')
xlabel('Stimulus Level s')
ylabel('Posterior Density')
xlim([-25,25])
%legend('C=0', 'C=1', 'Location', 'northwest')
%legend('boxoff')
title('(Proto-) Posterior')

subplot(2,3,3)
hold on
y = y1 .* l + y2 .* l;
y1f = y .* sum(y1 .* l);
y2f = y .* sum(y2 .* l);
patch(x, y2f, 'red', 'FaceAlpha', 0.3, 'LineWidth', 2, 'EdgeColor', 'red');
patch(x, y1f, 'blue', 'FaceAlpha', 0.3, 'LineWidth', 2, 'EdgeColor', 'blue');
[m1f, ~] = max(y1f);
[m2f, xf] = max(y2f);
xf = x(xf);
plot([-30,xf], [m1f, m1f], 'k:', 'LineWidth', 2);
plot([-30,xf], [m2f, m2f], 'k:', 'LineWidth', 2);
plot([xf,xf], [0, m1f], 'k:', 'LineWidth', 2);
plot([xf,xf], [0, m2f], 'k:', 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'Fontsize', 14)
set(gca, 'XTick', [0, x0, xf])
set(gca, 'YTick', [], 'TickDir', 'out')
set(gca, 'XTickLabel', {'0', '', '$\hat{s}$'}, 'TickLabelInterpreter','latex')
xlabel('Stimulus Level s')
ylabel('Factorized approximation q')
xlim([-25,25])
%legend('C=0', 'C=1', 'Location', 'northwest')
%legend('boxoff')
title('Factorized Approximation')


x0 = 3;
y1 = 2 * exp(-0.5 * x.^2 ./ sigma2.^2) ./ sqrt(2 *pi) ./ sigma2 .* (x <= 0);
y2 = 2 * exp(-0.5 * x.^2 ./ sigma2.^2) ./ sqrt(2 *pi) ./ sigma2 .* (x > 0);
l = exp(-0.5 * (x-x0).^2 ./ sigmaL.^2) ./ sqrt(2 *pi) ./ sigmaL;


subplot(2,3,4)
hold on
plot(x, y1, 'LineWidth', 2);
plot(x, y2, 'LineWidth', 2);
plot(x, l, '--k', 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'Fontsize', 14)
set(gca, 'XTick', [0, x0])
set(gca, 'YTick', [])
set(gca, 'XTickLabel', {'0', '$x$'}, 'TickLabelInterpreter','latex')
xlabel('Stimulus Level s')
ylabel('Density')
legend('$C=0$', '$C=1$', 'likelihood', 'Location', 'northwest', ...
    'Interpreter','latex')
legend('boxoff')
xlim([-25,25])

subplot(2,3,5)
hold on
patch(x, y2 .* l, 'red', 'FaceAlpha', 0.3, 'LineWidth', 2, 'EdgeColor', 'red');
patch(x, y1 .* l, 'blue', 'FaceAlpha', 0.3, 'LineWidth', 2, 'EdgeColor', 'blue');
[m1, x1] = max(y1 .* l);
[m2, x2] = max(y2 .* l);
x1 = x(x1);
x2 = x(x2);
plot([-30,x1], [m1, m1], 'k:', 'LineWidth', 2);
plot([-30,x2], [m2, m2], 'k:', 'LineWidth', 2);
plot([x1,x1], [0, m1], 'k:', 'LineWidth', 2);
plot([x2,x2], [0, m2], 'k:', 'LineWidth', 2);
set(gca, 'TickDir', 'out', 'Fontsize', 14)
set(gca, 'XTick', [x1, x2, x0])
set(gca, 'YTick', [])
set(gca, 'XTickLabel', {'$\hat{s}_0$', '$\hat{s}_1$', ''}, 'TickLabelInterpreter','latex')
xlabel('Stimulus Level s')
ylabel('Posterior Density')
xlim([-25,25])
%legend('C=0', 'C=1', 'Location', 'northwest')
%legend('boxoff')

% subplot(2,3,6)
% hold on
% y = y1 .* l + y2 .* l;
% y1f = y .* sum(y1 .* l);
% y2f = y .* sum(y2 .* l);
% patch(x, y2f, 'red', 'FaceAlpha', 0.3, 'LineWidth', 2, 'EdgeColor', 'red');
% patch(x, y1f, 'blue', 'FaceAlpha', 0.3, 'LineWidth', 2, 'EdgeColor', 'blue');
% [m1f, ~] = max(y1f);
% [m2f, xf] = max(y2f);
% xf = x(xf);
% plot([-30,xf], [m1f, m1f], 'k:', 'LineWidth', 2);
% plot([-30,xf], [m2f, m2f], 'k:', 'LineWidth', 2);
% plot([xf,xf], [0, m1f], 'k:', 'LineWidth', 2);
% plot([xf,xf], [0, m2f], 'k:', 'LineWidth', 2);
% set(gca, 'TickDir', 'out', 'Fontsize', 14)
% set(gca, 'XTick', [0, xf, x0])
% set(gca, 'YTick', [], 'TickDir', 'out')
% set(gca, 'XTickLabel', {'0', '$\hat{x}$', ''}, 'TickLabelInterpreter','latex')
% xlabel('Stimulus Level s')
% ylabel('Factorized approximation q')
xlim([-25,25])
%legend('C=0', 'C=1', 'Location', 'northwest')
%legend('boxoff')


set(gcf, 'Position', [250, 650, 1050, 500])

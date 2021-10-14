function plot_decision(type)

N = 100;

x1 = linspace(-50, 50, N);
x2 = linspace(-50, 50, N);

[x1_grid, x2_grid] = meshgrid(x1, x2);

x = cat(2,x1_grid(:), x2_grid(:));

sigmaNoise = 15;
sigma0 = 24;
lambda = 0;
beta0 = 0;
beta = 1;

if strcmp(type,'bayes')
    priorD = sigmaNoise .* sqrt(2*sigmaNoise.^2 + sigma0.^2) ./ (sigmaNoise.^2 + sigma0.^2);
    pd = log(priorD);

    d = pd - (x(:,1).^2 + x(:,2).^2)./2./(sigma0.^2 + sigmaNoise.^2);
    d = d  - (x(:,1) + x(:,2)).^2 .* sigma0.^2 ./ 2./ sigmaNoise.^2 ./ (2*sigma0.^2 + sigmaNoise.^2);
    d = d  + (x(:,1).^2 + x(:,2).^2)./2 ./ sigmaNoise.^2;
elseif strcmp(type,'PE')
    shatSeparate = (sigma0.^2 ./ (sigma0.^2 + sigmaNoise.^2)) .* x;
    shatSame = (sigma0.^2 ./ (2*sigma0.^2 + sigmaNoise.^2)) .* sum(x, 2);

    d = - 1 / 2 .* ((sum(shatSeparate.^2, 2) - shatSame.^2) ./ sigma0.^2 ...
        + (sum((x - shatSeparate).^2, 2) - sum((x - repmat(shatSame,[1,2])).^2, 2))./ sigmaNoise);
end
p = lambda./ 2 + (1-lambda)./(1+exp(beta0 + beta.*d));
p = reshape(p, N, N);

imagesc(x1,x2,p, [0,1])

%contour(x1, x2, p, linspace(0,1,11))
xlabel('Left Measurement [px]', 'FontSize', 24)
ylabel('Right Measurement [px]', 'FontSize', 24)
set(gca, 'FontSize', 18)
colorbar()
axis square
axis tight

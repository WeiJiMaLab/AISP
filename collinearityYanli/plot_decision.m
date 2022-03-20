function plot_decision(type, n_rep)

N = 100;
n_samples = 100;
if ~exist('n_rep', 'var') || isempty(n_rep)
    n_rep = 10;
end

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
    p = lambda./ 2 + (1-lambda)./(1+exp(beta0 + beta.*d));
elseif strcmp(type,'PE')
    shatSeparate = (sigma0.^2 ./ (sigma0.^2 + sigmaNoise.^2)) .* x;
    shatSame = (sigma0.^2 ./ (2*sigma0.^2 + sigmaNoise.^2)) .* sum(x, 2);

    d = - 1 / 2 .* ((sum(shatSeparate.^2, 2) - shatSame.^2) ./ sigma0.^2 ...
        + (sum((x - shatSeparate).^2, 2) - sum((x - repmat(shatSame,[1,2])).^2, 2))./ sigmaNoise);
    p = lambda./ 2 + (1-lambda)./(1+exp(beta0 + beta.*d));
elseif strcmp(type,'sample')
    p = zeros(size(x, 1), 1);
    for i = 1:n_rep
        evidence = zeros(size(x, 1), 2);
        s1_samp = sigma0 * repmat(randn(n_samples, size(x, 1), 1), [1,1,2]);
        s2_samp = sigma0 * randn(n_samples, size(x, 1), 2);
        evidence(:,1) = logsumexp(-sum((reshape(x,[1,size(x,1),2])-s1_samp).^2 ./ sigmaNoise.^2 ./2,3));
        evidence(:,2) = logsumexp(-sum((reshape(x,[1,size(x,1),2])-s2_samp).^2 ./ sigmaNoise.^2 ./2,3));
        d = evidence(:,2) - evidence(:,1);
        p = p +lambda./ 2 + (1-lambda)./(1+exp(beta0 + beta.*d));
    end
    p = p ./ n_rep;
end
p = reshape(p, N, N);

imagesc(x1,x2,p, [0,1])

%contour(x1, x2, p, linspace(0,1,11))
xlabel('Left Measurement [px]', 'FontSize', 24)
ylabel('Right Measurement [px]', 'FontSize', 24)
set(gca, 'FontSize', 18)
colorbar()
axis square
axis tight

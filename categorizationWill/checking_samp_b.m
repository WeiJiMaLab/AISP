

xx = linspace(-20,20,250);
sigmaNoisex = 4;

d = 1/2*log((sigma2.^2+sigmaNoisex.^2)./(sigma1.^2+sigmaNoisex.^2))- ...
    xx.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2+sigmaNoisex.^2)./(sigma2^2+sigmaNoisex^2);

d_samp = zeros(4,250);
k = 0;
for n_samples = [10,100,1000,10000]
    k = k + 1;
    p = 0;
    for i = 1:n_rep
        s1_samp = sigma1 * randn(size(xx, 1), size(xx, 2), n_samples);
        s2_samp = sigma2 * randn(size(xx, 1), size(xx, 2), n_samples);
        evidence1 = logsumexp(- (s1_samp-repmat(xx, 1, 1, n_samples)).^2 ./ 2 ./ repmat(sigmaNoisex .^2, 1, 1, n_samples), 3);
        evidence2 = logsumexp(- (s2_samp-repmat(xx, 1, 1, n_samples)).^2 ./ 2 ./ repmat(sigmaNoisex .^2, 1, 1, n_samples), 3);
        d_samp(k,:) = d_samp(k,:) + evidence1 - evidence2;
    end
end
d_samp = d_samp ./ n_rep;

x_int = linspace(-50,50,2000);
y1 = exp(-x_int.^2 ./ 2 ./ sigma1.^2) / sqrt(2*pi) / sigma1;
y2 = exp(-x_int.^2 ./ 2 ./ sigma2.^2) / sqrt(2*pi) / sigma2;

% clf 
% plot(x_int, y1)
% hold on
% plot(x_int, y2)

dist = (x_int - xx') .^2;
ev1 = sum(y1 .* exp(-dist./2./(sigmaNoisex.^2)), 2);
ev2 = sum(y2 .* exp(-dist./2./(sigmaNoisex.^2)), 2);

clf
subplot(1,2,1);
plot(xx, d);
hold on
plot(xx,log(ev1./ev2));
plot(xx, d_samp);
subplot(1,2,2);
plot(xx, 1./(1+exp(-d)));
hold on
plot(xx,1./(1+exp(-log(ev1./ev2))));
plot(xx, 1./(1+exp(-d_samp)));
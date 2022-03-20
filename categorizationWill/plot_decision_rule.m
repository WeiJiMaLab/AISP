function plot_decision_rule(type, n_rep, n_samples)

if ~exist('n_rep','var') || isempty(n_rep)
    n_rep = 10;
end

x = linspace(-20,20,250);
sigmaNoise = exp(linspace(log(0.099),log(20),251));

[xx,sigmaNoisex] = meshgrid(x,sigmaNoise);

sigma1 = 3;
sigma2 = 12;

switch type
    case 'Freq'
        d = 1/2*log((sigma2.^2)./(sigma1.^2))- ...
            xx.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2+sigmaNoisex.^2)./(sigma2.^2+sigmaNoisex.^2);
        p = 1./(1+exp(-d));
    case 'Freq2'
        d = 1/2*log((sigma2.^2)./(sigma1.^2))- ...
            xx.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2)./(sigma2^2);
        p = 1./(1+exp(-d));
    case 'Freq3'
        w1 = normpdf(xx,0,sigma1^2+sigmaNoisex);
        w2 = normpdf(xx,0,sigma2^2+sigmaNoisex);
        w1n = w1./(w1+w2);
        w2n = w2./(w1+w2);
        s = w1n .* sigma1./(sigma1+sigmaNoisex).*xx+w2n .* sigma2./(sigma2+sigmaNoisex).*xx;
        d = 1/2*log((sigma2.^2)./(sigma1.^2))- ...
            s.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2)./(sigma2^2);
        p = 1./(1+exp(-d));
    case 'Bayes'
        d = 1/2*log((sigma2.^2+sigmaNoisex.^2)./(sigma1.^2+sigmaNoisex.^2))- ...
            xx.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2+sigmaNoisex.^2)./(sigma2.^2+sigmaNoisex.^2);
        p = 1./(1+exp(-d));
    case 'Var'
        tol = 0.00001;
        q1 = zeros(size(xx));
        q1new = 0.5 * ones(size(xx));
        while any(abs(q1new-q1)>tol)
            q1 = q1new;
            sigmanew = (sigma1.^2.*sigma2.^2)./(q1*sigma1.^2.*sigmaNoisex.^2 + ...
                (1-q1).*sigma2.^2.*sigmaNoisex.^2+sigma1.^2*sigma2.^2);
            shat = sigmanew.*xx;
            svar = sigmanew.*sigmaNoisex.^2;
            z = shat.^2+svar;
            q1new = exp(-z./2./sigma1.^2)./sigma1./(exp(-z./2./sigma1.^2)./sigma1 + exp(-z./2./sigma2.^2)./sigma2);
        end
        q1 = q1new;
        d = log(q1)-log(1-q1);
        p = 1./(1+exp(-d));
    case 'Sample'
        if ~exist('n_samples','var') || isempty(n_samples)
            n_samples = 60;% mean over all fits was 64.7716
        end
        p = 0;
        for i = 1:n_rep
            s1_samp = sigma1 * randn(size(xx, 1), size(xx, 2), n_samples);
            s2_samp = sigma2 * randn(size(xx, 1), size(xx, 2), n_samples);
            evidence1 = logsumexp(- (s1_samp-repmat(xx, 1, 1, n_samples)).^2 ./ 2 ./ repmat(sigmaNoisex .^2, 1, 1, n_samples), 3);
            evidence2 = logsumexp(- (s2_samp-repmat(xx, 1, 1, n_samples)).^2 ./ 2 ./ repmat(sigmaNoisex .^2, 1, 1, n_samples), 3);
            d = evidence1 - evidence2;
            p = p + 1./(1+exp(-d));
        end
        p = p / n_rep;
    case 'cssample'
        if ~exist('n_samples','var') || isempty(n_samples)
            n_samples = 150; % mean over all fits was 154.8879
        end
        sigma12 = [sigma1, sigma2];
        sigmaNoise = sigmaNoisex(:);
        lsn = log(sigmaNoise(:));
        x = xx(:);
        p = 0;
        n_trials = numel(xx);
        for i = 1:n_rep
            C_samp = randi(2, n_trials, 1);
            s_samp = sigma12(C_samp)' .* randn(n_trials, 1);
            s_s = nan(n_trials, n_samples);
            C = nan(n_trials, n_samples);
            
            l_samp = - (s_samp-x).^2 ./ 2 ./ sigmaNoise .^2 - lsn;
            % sampling
            for j = 1:n_samples
                % sample C
                C_new = randi(2, n_trials, 1);
                % sample s
                s_samp_new = sigma12(C_new)' .* randn(n_trials, 1);
                l_samp_new = - (s_samp_new-x).^2 ./ 2 ./ sigmaNoise .^2 - lsn;
                accept = rand(n_trials,1) < exp(l_samp_new-l_samp);
                s_samp(accept) = s_samp_new(accept);
                l_samp(accept) = l_samp_new(accept);
                C_samp(accept) = C_new(accept);
                C(:,j) = C_samp;
                s_s(:,j) = s_samp;
            end
            % 1 is the wide category = 2 here
            d = log(sum(C==1, 2) + 1) - log(sum(C==2, 2) + 1);
            p = p + 1./(1+exp(-d));
        end
        p = p / n_rep;
        p = reshape(p, size(xx));
end

imagesc(x,sigmaNoise,p,[0,1])
set(gca,'YScale','log')
ylim([min(sigmaNoise),max(sigmaNoise)])
colormap(parula)
hc = colorbar();
set(hc, 'Ticks', [0, 0.5, 1], 'Color', 'k')
xlabel('Measurement x','FontSize',16, 'Color', 'k')
ylabel('\sigma_n Noise standard deviation', 'FontSize',16, 'Color', 'k')
xticks([-20,-10,0,10,20])
yticks([0.1,1,10])
yticklabels({'0.1', '1', '10'})
set(gca, 'TickDir', 'out', 'FontSize', 18, 'LineWidth', 2, 'xcolor', 'k', 'ycolor','k')
box off
set(gcf, 'Position', [560,750,560,250])

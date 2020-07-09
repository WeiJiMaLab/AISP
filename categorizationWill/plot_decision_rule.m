function plot_decision_rule(type)

x = linspace(-20,20,250);
sigmaNoise = exp(linspace(log(0.1),log(20),250));

[xx,sigmaNoisex] = meshgrid(x,sigmaNoise);

sigma1 = 3;
sigma2 = 12;

switch type
    case 'Freq'
        d = 1/2*log((sigma2.^2)./(sigma1.^2))- ...
            xx.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2+sigmaNoisex.^2)./(sigma2^2+sigmaNoisex^2);
    case 'Freq2'
        d = 1/2*log((sigma2.^2)./(sigma1.^2))- ...
            xx.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2)./(sigma2^2);
    case 'Freq3'
        w1 = normpdf(xx,0,sigma1^2+sigmaNoisex);
        w2 = normpdf(xx,0,sigma2^2+sigmaNoisex);
        w1n = w1./(w1+w2);
        w2n = w2./(w1+w2);
        s = w1n .* sigma1./(sigma1+sigmaNoisex).*xx+w2n .* sigma2./(sigma2+sigmaNoisex).*xx;   
        d = 1/2*log((sigma2.^2)./(sigma1.^2))- ...
            s.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2)./(sigma2^2);
    case 'Bayes'
        d = 1/2*log((sigma2.^2+sigmaNoisex.^2)./(sigma1.^2+sigmaNoisex.^2))- ...
            xx.^2/2 * (sigma2.^2 - sigma1.^2)./(sigma1.^2+sigmaNoisex.^2)./(sigma2^2+sigmaNoisex^2);
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
    
end

imagesc(x,sigmaNoise,1./(1+exp(-d)),[0,1])
set(gca,'YScale','log')
colormap(parula)
colorbar()
xlabel('Measurement x','FontSize',14)
ylabel('\sigma_n Noise standard deviation','FontSize',14)
xticks([])
set(gca, 'TickDir', 'out')
box off
set(gcf, 'Position', [560,750,560,190])

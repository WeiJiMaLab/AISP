function plot_decision_rule(type)

x = linspace(-25,25,250);
sigmaNoise = linspace(0,10,250);

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
    
end

imagesc(x,sigmaNoise,1./(1+exp(-d)),[0,1])
colormap(gray)
colorbar()
xlabel('Measurement x','FontSize',14)
ylabel('\sigma_n Noise standard deviation','FontSize',14)
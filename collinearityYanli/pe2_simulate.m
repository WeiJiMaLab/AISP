function responses = pe2_simulate(stim, pars)
% function responses = pe2_simulate(data,sigmas,beta,lambda)
% simulates respones of the bayesian observer.
% stim should be a matrix of size trials x 3 with columns:
%     mean 
%     stimulus s1
%     stimulus s2 
% sigmas should be a 3 element vector of noise standard deviations
% beta should be the constant and slope for the logistic link
% lambda is the lapserate

persistent optimal_criterion
persistent saved_pars

sigmas = pars(1:4);
beta0 = pars(5);
beta = pars(6);
lambda = pars(7);
sigma0 = 24;

Nsim = 500000;
if isempty(saved_pars) || any(saved_pars ~=pars)
    optimal_criterion = zeros(4);
    for i_sigma = 1:4
        x_sim0 = sqrt(sigmas(i_sigma).^2 + sigma0.^2) * randn(Nsim,2);
        shatSeparate0 = (sigma0.^2 ./ (sigma0.^2 + sigmas(i_sigma).^2)) .* x_sim0;
        shatSame0 = (sigma0.^2 ./ (2*sigma0.^2 + sigmas(i_sigma).^2)) .* sum(x_sim0, 2);
        d0 = - 1 / 2 .* ((sum(shatSeparate0.^2, 2) - shatSame0.^2) ./ sigma0.^2 ...
            + (sum((x_sim0 - shatSeparate0).^2, 2) - sum((x_sim0 - repmat(shatSame0,[1,2])).^2, 2))./ sigmas(i_sigma));
        x_sim1 = repmat(sigma0 * randn(Nsim,1), [1,2]);
        x_sim1 = x_sim1 + sigmas(i_sigma) * randn(Nsim,2);
        shatSeparate1 = (sigma0.^2 ./ (sigma0.^2 + sigmas(i_sigma).^2)) .* x_sim1;
        shatSame1 = (sigma0.^2 ./ (2*sigma0.^2 + sigmas(i_sigma).^2)) .* sum(x_sim1, 2);
        d1 = - 1 / 2 .* ((sum(shatSeparate1.^2, 2) - shatSame1.^2) ./ sigma0.^2 ...
            + (sum((x_sim1 - shatSeparate1).^2, 2) - sum((x_sim1 - repmat(shatSame1,[1,2])).^2, 2))./ sigmas(i_sigma));
        
        % Bisection method homebrew:
        crit0 = -20;
        crit2 = 100;
        crit1 = crit0 + (crit2-crit0) / ( 3 + sqrt(5) ) * 2;
        f0 = sum(d1 < crit0) + sum(d0 > crit0);
        f1 = sum(d1 < crit1) + sum(d0 > crit1);
        f2 = sum(d1 < crit2) + sum(d0 > crit2);
        while (crit2-crit0) > 0.001
            if abs(crit2-crit1) > abs(crit1-crit0)
                crit_new = crit1 + (crit2-crit1) / ( 3 + sqrt(5) ) * 2;
                f_new = sum(d1 < crit_new) + sum(d0 > crit_new);
                if f_new < f1
                    f2 = f_new;
                    crit2 = crit_new;
                else
                    f0 = f1;
                    crit0 = crit1;
                    f1 = f_new;
                    crit1 = crit_new;
                end
            else
                crit_new = crit1 + (crit0-crit1) / ( 3 + sqrt(5) ) * 2;
                f_new = sum(d1 < crit_new) + sum(d0 > crit_new);
                if f_new < f1
                    f0 = f_new;
                    crit0 = crit_new;
                else
                    f2 = f1;
                    crit2 = crit1;
                    f1 = f_new;
                    crit1 = crit_new;
                end
            end
        end
        optimal_criterion(i_sigma) = crit1;
    end
    saved_pars = pars;
end


s = stim(:, [2,3]) - stim(:, 1);

sigmaNoise = zeros(size(stim, 1), 1);
sigmaNoise(stim(:,1)==0) = sigmas(1);
sigmaNoise(stim(:,1)==240) = sigmas(2);
sigmaNoise(stim(:,1)==480) = sigmas(3);
sigmaNoise(stim(:,1)==840) = sigmas(4);
oc = zeros(size(stim, 1), 1);
oc(stim(:,1)==0) = optimal_criterion(1);
oc(stim(:,1)==240) = optimal_criterion(2);
oc(stim(:,1)==480) = optimal_criterion(3);
oc(stim(:,1)==840) = optimal_criterion(4);
x = s + repmat(sigmaNoise, [1,2]) .* randn(size(stim, 1), 2);

shatSeparate = (sigma0.^2 ./ (sigma0.^2 + sigmaNoise.^2)) .* x;
shatSame = (sigma0.^2 ./ (2*sigma0.^2 + sigmaNoise.^2)) .* sum(x, 2);

d = - 1 / 2 .* ((sum(shatSeparate.^2, 2) - shatSame.^2) ./ sigma0.^2 ...
    + (sum((x - shatSeparate).^2, 2) - sum((x - repmat(shatSame,[1,2])).^2, 2))./ sigmaNoise) ...
    - oc; % subtract optimal criterion
p = lambda./2 + (1-lambda)./(1+exp(beta0 + beta*d));
responses = (rand(size(stim, 1), 1) < p);

function responses = sample_simulate(stim, pars)
% function responses = bayes_simulate(data,pars)
% simulates respones of the bayesian observer.
% stim should be a matrix of size trials x 3 with columns:
%     mean 
%     stimulus s1
%     stimulus s2

sigmas = pars(1:4);
beta0 = pars(5);
beta = pars(6);
lambda = pars(7);
if length(pars) > 7
    n_samples = floor(pars(8));
    if n_samples > 10000
        n_samples = 10000;
        warning('n_samples reduced to 10000 to avoid extreme computation')
    elseif rand > mod(n_samples, 1)
        n_samples = floor(n_samples);
    else
        n_samples = ceil(n_samples);
    end
else
    n_samples = 1;
end
sigma0 = 24;

s = stim(:, [2,3]) - stim(:, 1);

sigmaNoise = zeros(size(stim, 1), 1);
sigmaNoise(stim(:,1)==0) = sigmas(1);
sigmaNoise(stim(:,1)==240) = sigmas(2);
sigmaNoise(stim(:,1)==480) = sigmas(3);
sigmaNoise(stim(:,1)==840) = sigmas(4);

x = s + repmat(sigmaNoise, [1,2]) .* randn(size(stim, 1), 2);
evidence = zeros(size(stim, 1), 2);
s1_samp = sigma0 * repmat(randn(n_samples, size(stim, 1), 1), [1,1,2]);
s2_samp = sigma0 * randn(n_samples, size(stim, 1), 2);
evidence(:,1) = logsumexp(-sum((reshape(x,[1,size(x,1),2])-s1_samp).^2 ./ reshape(sigmaNoise, [1,size(x,1),1]).^2 ./2,3));
evidence(:,2) = logsumexp(-sum((reshape(x,[1,size(x,1),2])-s2_samp).^2 ./ reshape(sigmaNoise, [1,size(x,1),1]).^2 ./2,3));
d = evidence(:,2) - evidence(:,1);
p = 1 ./(1 + exp(beta0 + beta.*d));
p(evidence(:,1) == 0) = 0;
p(evidence(:,2) == 0) = 1;
p = lambda./ 2 + (1-lambda) .* p;
responses = (rand(size(stim, 1), 1) < p);

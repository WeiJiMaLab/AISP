function test_sample_bayes(n_samples)

data = getData;
pars = load('pars/pars_Bayes_1_1.mat').pars;

stim = data(:, [2, 5:6]);


sigmas = pars(1:4);
beta0 = pars(5);
beta = pars(6);
lambda = pars(7);
sigma0 = 24;

priorD = sigmas .* sqrt(sigmas.^2 + 2*sigma0.^2) ./ (sigmas.^2 + sigma0.^2);
priorD = log(priorD);
        
s = stim(:, [2,3]) - stim(:, 1);

sigmaNoise = zeros(size(stim, 1), 1);
sigmaNoise(stim(:,1)==0) = sigmas(1);
sigmaNoise(stim(:,1)==240) = sigmas(2);
sigmaNoise(stim(:,1)==480) = sigmas(3);
sigmaNoise(stim(:,1)==840) = sigmas(4);
pd = zeros(size(stim, 1), 1);
pd(stim(:,1)==0) = priorD(1);
pd(stim(:,1)==240) = priorD(2);
pd(stim(:,1)==480) = priorD(3);
pd(stim(:,1)==840) = priorD(4);

x = s + repmat(sigmaNoise, [1,2]) .* randn(size(stim, 1), 2);
% old wrong formulas:
%d = pd - sigma0.^2 ./ 2 ./ sigmaNoise.^2 ./ (2 .* sigmaNoise.^2 + sigma0.^2) .* x(:,1) .* x(:,2);
%d = d + sigma0.^2 ./ (sigma0.^2 + sigmaNoise.^2) ./ (4 .* sigmaNoise.^2 + 2 * sigma0.^2) .* (x(:,1).^2 + x(:,2).^2);
d = pd - (x(:,1).^2 + x(:,2).^2)./2./(sigma0.^2 + sigmaNoise.^2);
d = d  - (x(:,1) + x(:,2)).^2 .* sigma0.^2 ./ 2./ sigmaNoise.^2 ./ (2*sigma0.^2 + sigmaNoise.^2);
d_bayes = d  + (x(:,1).^2 + x(:,2).^2)./2 ./ sigmaNoise.^2;



evidence = zeros(size(stim, 1), 2);
s1_samp = sigma0 * repmat(randn(n_samples, size(stim, 1), 1), [1,1,2]);
s2_samp = sigma0 * randn(n_samples, size(stim, 1), 2);
evidence(:,1) = logsumexp(-sum((reshape(x,[1,size(x,1),2])-s1_samp).^2 ./ reshape(sigmaNoise, [1,size(x,1),1]).^2 ./2,3));
evidence(:,2) = logsumexp(-sum((reshape(x,[1,size(x,1),2])-s2_samp).^2 ./ reshape(sigmaNoise, [1,size(x,1),1]).^2 ./2,3));
%evidence(:, 1) = evidence(:,1) + exp(sum(-(x-s1_samp).^2 ./ sigmaNoise.^2 ./ 2, 2));
%evidence(:, 2) = evidence(:,2) + exp(sum(-(x-s2_samp).^2 ./ sigmaNoise.^2 ./ 2, 2));
% evidence = evidence ./ n_samples; % unnecessary due to ratio in next step
d_samp = evidence(:,2) - evidence(:,1);
d_samp(d_samp>2000) = 2000;
d_samp(d_samp<-2000) = -2000;
d_bayes(d_bayes>2000) = 2000;
d_bayes(d_bayes<-2000) = -2000;
subplot(1,2,1)
hist(d_samp-d_bayes)
subplot(1,2,2)
plot(d_samp, d_bayes, 'k.')

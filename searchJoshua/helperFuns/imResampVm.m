function samples = imResampVm(mu, kappa, nSamples)
% Use sampling importance resampling, to approximately sample from the 
% von Mises distribution

% INPUT
% mu: scalar. Centre of the distribution
% kappa: scalar. Concentration parameter
% nSamples: scalar. How many samples to draw?

% JCT, 2021

persistent evenSpaced

assert(length(mu) == 1)
assert(length(kappa) == 1)

if isempty(evenSpaced)
    evenSpaced = -pi : 0.001 : pi;
end
assocProbs = circ_vmpdf(evenSpaced, mu, kappa);

samples = randsample(evenSpaced, nSamples, true, assocProbs);



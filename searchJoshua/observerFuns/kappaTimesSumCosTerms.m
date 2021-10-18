function res = kappaTimesSumCosTerms(percept, stimSamples, kappa_x)
% A calculation frequently encountered in the derivations

summed = sum(cos(percept - stimSamples), 2, 'omitnan');
res = kappa_x .* summed;
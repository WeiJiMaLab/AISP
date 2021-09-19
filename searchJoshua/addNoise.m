function noisy = addNoise(noiseless, kappas)
% Add VM noise to the values in the 2D matrix noiseless. Values in
% noiseless that are nan remain as nan. kappas is a vector as long as the
% first dimension of noiseless, and determines the kappa of the noise used
% for each row. Values are mapped back into the range -2pi to pi

% JCT, 2021

assert(ndims(noiseless)==2)
assert(ndims(kappas)==2)
assert(size(kappas, 2)==1)
assert(size(kappas, 1)==size(noiseless, 1))

inactive = isnan(noiseless);
kappaShaped = repmat(kappas, 1, size(noiseless, 2));
kappaShaped(inactive) = nan;

addedNoise = qrandvm(0, kappaShaped, size(noiseless));
noisy = addedNoise + noiseless;
noisy = vS_mapBackInRange(noisy, -pi, pi);

assert(isequal(isnan(noisy), isnan(noiseless)))
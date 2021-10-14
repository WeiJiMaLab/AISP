function noisy = addNoise(noiseless, kappas)
% Add VM noise to the values in the matrix noiseless. Values in
% noiseless that are nan remain as nan. Values are mapped back into the 
% range -2pi to pi

% INPUT
% noiseless: matrix (may be 3D)
% kappas: vector as long as the first dimension of noiseless, and 
%   determines the kappa of the noise used for each row. 

% JCT, 2021

assert(ndims(kappas)==2)
assert(size(kappas, 2)==1)
assert(size(kappas, 1)==size(noiseless, 1))

nlShape = size(noiseless);
repFactors = [1, nlShape(2:end)];
kappaShaped = repmat(kappas, repFactors);
assert(isequal(size(kappaShaped), nlShape))

inactive = isnan(noiseless);
kappaShaped(inactive) = nan;

addedNoise = qrandvm(0, kappaShaped, nlShape);
noisy = addedNoise + noiseless;
noisy = vS_mapBackInRange(noisy, -pi, pi);

assert(isequal(isnan(noisy), isnan(noiseless)))
function dvFunsInputChecks(percept, mu_s, nItems, kappa_s, kappa_x)

assert(mu_s == 0)
assert(not(any(all(isnan(percept),2)))) % We use summation with the 
% 'omitnan' flag, and this would incorrectly give an numeric result 
% for a row that was all nans

if size(percept, 2) > 8; error('Bug'); end
assert(length(size(percept)) == 2)

inputVectors = {nItems, kappa_s, kappa_x};

for iInputVec = 1 : length(inputVectors)
    vecSize = size(inputVectors{iInputVec});
    if (length(vecSize) ~= 2) || (vecSize(2) ~= 1)
        error('Bug')
    end
end
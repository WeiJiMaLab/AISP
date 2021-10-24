function thisSize = find3Dsize(thisMatrix)
% Check thisMatrix is 3D or less, then return size as a three element
% vector. Differs from size(thisMatrix) which would return a two element
% vector if thisMatrix was 2- or 1D.

assert(length(size(thisMatrix)) <= 3)
thisSize = [size(thisMatrix, 1), size(thisMatrix, 2), size(thisMatrix, 3)];

end
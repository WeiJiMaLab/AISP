function tidyMat = sendNansToEnd(mat)
% In each row, indedpendently, move all the nans to the end of the row

% INPUT
% mat: 2D matrix

% OUTPUT
% tidyMat: 2D matrix

% JCT, 2021, setup

assert(length(size(mat)) == 2)

tidyMat = nan(size(mat));
activeLocs = ~isnan(mat);
numActiv = sum(activeLocs, 2);
for iRow = 1 : size(mat, 1)
    thisRow = mat(iRow, :);
    tidyMat(iRow, 1:numActiv(iRow)) = thisRow(activeLocs(iRow, :));
end

assert(isequal(size(mat), size(tidyMat)))


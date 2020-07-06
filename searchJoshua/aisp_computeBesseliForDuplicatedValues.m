function besseliVals = aisp_computeBesseliForDuplicatedValues(vals)
% Compute besseli of order zero. Useful when vals contains many duplicated
% values. In this case compuation may be sped up by this function. 

uniqueVals = unique(vals);
uniqueResults = log(besseli(0, uniqueVals));

besseliVals = NaN(size(vals));

for iVal = 1 : length(uniqueVals)    
    besseliVals(vals == uniqueVals(iVal)) = uniqueResults(iVal);
    
end
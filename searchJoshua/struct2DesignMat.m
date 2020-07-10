function out = struct2DesignMat(in, direction)
% Converts a structure continaing the trial details to a design matrix 
%(direction=='to matrix'), or converts back to a structure (direction=='to struct')

relevantFields = {'SetSize', 'Orientation', 'KappaS', 'SetSizeCond'};
relevantCols = {1, 2:7, 8, 9};

if direction  == 'to matrix'
    out = nan(size(DatSubj.Response, 1), relevantCols{end}(end));
    
    for iF = 1 : length(relevantFields)
        dataToAdd = in.(relevantFields{iF});
        out(:, relevantCols{iF}) = dataToAdd;
    end
    
elseif direction  == 'to struct'
    for iF = 1 : length(relevantFields)
        out.(relevantFields{iF}) = in(:, relevantCols{iF});
    end
    
else
    error('Incorrect use of inputs')
end
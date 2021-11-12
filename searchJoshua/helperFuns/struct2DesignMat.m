function out = struct2DesignMat(in, direction, sortNans)
% Converts a structure continaing the trial details to a design matrix 
%(direction=='to matrix'), or converts back to a structure (direction=='to struct')

% INPUT
% sortNans: bool. If true, when running in direction 'to matrix', for 
%  any fields in 'in' that contain 2D matricies, for each row of these
%  matricies independently, send the nans to the end of the row.

relevantFields = {'SetSize', 'Orientation', 'KappaS', 'SetSizeCond'};
relevantCols = {1, 2:7, 8, 9};

if strcmp(direction, 'to matrix')
    out = nan(size(in.Response, 1), relevantCols{end}(end));
    
    for iF = 1 : length(relevantFields)
        dataToAdd = in.(relevantFields{iF});
        
        if sortNans && (length(size(dataToAdd)) == 2)
            dataToAdd = sendNansToEnd(dataToAdd);
        end
        
        out(:, relevantCols{iF}) = dataToAdd;
    end
    
elseif strcmp(direction, 'to struct')
    for iF = 1 : length(relevantFields)
        out.(relevantFields{iF}) = in(:, relevantCols{iF});
    end
    
else
    error('Incorrect use of inputs')
end
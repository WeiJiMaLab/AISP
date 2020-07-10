function out = paramVec2Struct(in, direction)
% Converts a vector of the parameters into a structure (direction=='to struct')

relevantFields = {'Kappa_x', 'Beta0', 'Beta1', 'LapseRate'};
relevantCols = {1:4, 5, 6, 7};

if direction == 'to struct'
    if (length(size(in))~=2) || (size(in, 1)~=1); error('Bug'); end
    
    for iF = 1 : length(relevantFields)
        out.(relevantFields{iF}) = in(relevantCols{iF});
    end
else
    error('Inocrrect use of inputs')
end
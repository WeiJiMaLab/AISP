function out = paramVec2Struct(in, direction)
% Converts a vector of the parameters into a structure (direction=='to struct')

[relevantFields, relevantCols] = findParamOrder();

if direction == 'to struct'
    if (length(size(in))~=2) || (size(in, 1)~=1); error('Bug'); end
    
    for iF = 1 : length(relevantFields)
        out.(relevantFields{iF}) = in(relevantCols{iF})';
    end
else
    error('Inocrrect use of inputs')
end
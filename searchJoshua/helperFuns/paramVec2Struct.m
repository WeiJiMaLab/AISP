function out = paramVec2Struct(in, modelName, direction)
% Converts a vector of the parameters into a structure (direction=='to struct')

if any(strcmp(modelName, {'impSamp', 'jointPostSamp'}))
    incSamplesParam = true;
else
    assert(any(strcmp(modelName, {'Bayes', 'PE', 'PE2', 'PE_imagineL'})))
    incSamplesParam = false;
end

[relevantFields, relevantCols] = findParamOrder(incSamplesParam);

if strcmp(direction, 'to struct')
    if (length(size(in))~=2) || (size(in, 1)~=1); error('Bug'); end
    
    used = false(size(in));
    for iF = 1 : length(relevantFields)
        out.(relevantFields{iF}) = in(relevantCols{iF})';
        
        assert(~any(used(relevantCols{iF})))
        used(relevantCols{iF}) = true;
    end
    assert(all(used(:)))
else
    error('Inocrrect use of inputs')
end
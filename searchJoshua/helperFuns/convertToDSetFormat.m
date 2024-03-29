function DSet = convertToDSetFormat(DSet, parsDir, Config)
% Takes the dataset DSet, and the fitting results in parsDir, and combines them
% to produce a single DSet which contains key fitting results in the format used
% by the modellingTools repository

% INPUT
% Config: struct. Has the following fields...
%   ModelLabel: Cell array of labels to use for each model
%   ModelList: Cell array of names of the models, as they were named during
%       fitting
%   Nreps: The minimim number of repitions for each model that was 
%       used during fitting.

if isfield(DSet.P, 'Models')
    error('DSet should not yet include fitting results.')
end

for iModel = 1 : length(Config.ModelList)
    for iPtpnt = 1 : length(DSet.P)
        DSet.P(iPtpnt).Models(iModel).Settings.ModelName = Config.ModelList{iModel};
    end
    
    files = findModelAssocFitFiles(parsDir, Config.ModelList{iModel});
    numFits = zeros(length(DSet.P), 1);
    
    for iFile = 1:length(files)
        f = load(fullfile(files(iFile).folder,files(iFile).name));
        fparts = split(files(iFile).name,{'_','.'});
        iPtpnt = str2double(fparts{end-2});
        assert(iPtpnt <= length(DSet.P))
        numFits(iPtpnt) = numFits(iPtpnt) +1;
        
        DSet.P(iPtpnt).Models(iModel).Fits(numFits(iPtpnt)).LL = -f.nLogL;
        DSet.P(iPtpnt).Models(iModel).Fits(numFits(iPtpnt)).Params ...
            = paramVec2Struct(f.pars, Config.ModelList{iModel}, ...
                                'to struct');
    end
end
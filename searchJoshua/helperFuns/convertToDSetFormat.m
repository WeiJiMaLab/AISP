function DSet = convertToDSetFormat(DSet, parsDir, configFile)
% Takes the dataset DSet, and the fitting results in parsDir, and combines them
% to produce a single DSet which contains key fitting results in the format used
% by the modellingTools repository

% INPUT
% configFile: string. File path to matlab file to be loaded.

warning('Check this function')

if isfield(DSet.P, 'Models')
    error('DSet should not yet include fitting results.')
end


Config = load(configFile);

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
        numFits(iPtpnt) = numFits(iPtpnt) +1;
        
        DSet.P(iPtpnt).Models(iModel).Fits(numFits(iPtpnt)).LL = -f.nLogL;
        error('Inputs to paramVec2Struct changed')
        DSet.P(iPtpnt).Models(iModel).Fits(numFits(iPtpnt)).Params ...
            = paramVec2Struct(f.pars, 'to struct');
    end
end
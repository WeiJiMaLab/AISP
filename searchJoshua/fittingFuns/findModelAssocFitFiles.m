function files = findModelAssocFitFiles(parsDir, modelName)
% Look in the directory parsDir, and find all the fit results files
% associated with modelName

% INPUT
% parsDir: String.
% modelName: String. Entry in Config.ModelList.

% OUTPUT
% files: Struct. Same format as produced by MATLAB function 'dir'.


fname = [parsDir '\pars_*_' modelName];
files = dir([fname,'_*']);


% With all the candidate files, run a stricter test to ensure we have
% located only those files we really want. This became necessary when I
% unfortunately started using a model called 'PE_imagineL', with the
% underscore messing the previous system up.
remove = zeros(length(files), 1);
for iFile = 1:length(files)
    
    pattern = [modelName '_\d*_\d*.mat'];
    match = regexp(files(iFile).name, pattern, 'once');
    if isempty(match); remove(iFile)=1; end
end
files(logical(remove))=[];
disp(['Model ' modelName ...
    ' files found: ' num2str(length(files))])
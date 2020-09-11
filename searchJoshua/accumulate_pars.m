function accumulate_pars(dataDir, parsDir, firstOrAll)

% INPUT
% dataDir: Directory containing the original unfitted dataset
% parsDir: Where to find the fitted params, and where to save the files with all
% the parameters collected
% firstOrAll: 'first' or 'all'. Collect the parameters associated with the first
% round of fitting only (up to Config.Nreps), or collect all parameters
% including any that were later scheduled using cluster_fcn_fancy

[~, Nptpnts] = getData(dataDir);
Config = load('Config.mat');
Nreps = Config.Nreps;

for itype = 1 : length(Config.ModelList)
    fname = [parsDir '\pars_*_' Config.ModelList{itype}];
    shortName = [parsDir '\pars_' Config.ModelList{itype}];
    files = dir([fname,'_*']);
    
    % Load the first file just to find out the number of parameters
    f = load(fullfile(files(1).folder,files(1).name));
    Npars = length(f.pars);
    pars = nan(Nptpnts,Npars,Nreps);
    nLogLs = nan(Nptpnts,Nreps);
    
    for iFile = 1:length(files)
        f = load(fullfile(files(iFile).folder,files(iFile).name));
        fparts = split(files(iFile).name,{'_','.'});
        iPtpnt = str2double(fparts{end-2});
        iRep = str2double(fparts{end-1});
        if strcmp(firstOrAll, 'first')
            if iRep <= Config.Nreps
                if (~isnan(nLogLs(iPtpnt,iRep))) ...
                        || (any(~isnan(pars(iPtpnt, :, iRep))))
                    error('There already seems to be data here.')
                end
                nLogLs(iPtpnt,iRep) = f.nLogL;
                pars(iPtpnt,:,iRep) = f.pars;
            end
        else
            error('Not coded up yet')
        end
    end
    
    save([shortName,'.mat'],'pars','nLogLs')
end

assert(~any(isnan(pars(:))))
assert(~any(isnan(nLogLs(:))))


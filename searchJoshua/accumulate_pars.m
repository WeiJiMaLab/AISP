function accumulate_pars(dataDir, parsDir, firstOrAll, varargin)

% INPUT
% dataDir: Directory containing the original unfitted dataset
% parsDir: Where to find the fitted params, and where to save the files with all
% the parameters collected
% firstOrAll: 'first' or 'all'. Collect the parameters associated with the first
% round of fitting only (up to Config.Nreps), or collect all parameters
% including any that were later scheduled using cluster_fcn_fancy
% varargin: Boolean. Override error if there appear to be duplicate fits? Use with
% caution, probably indicates a bug.

if ~isempty(varargin)
    overrideDupError = varargin{1};
else
    overrideDupError = false;
end

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
        
        if strcmp(firstOrAll, 'first') && (iRep > Config.Nreps)
            continue % We only want to process the first fits
        end
        
        % Do we need to expand the pars and nLogLs arrays?
        if strcmp(firstOrAll, 'all')
            [pars, nLogLs] = expandArrays(pars, nLogLs, iRep, ...
                Nptpnts, Npars);
        end
            
        % Add in found values
        if (~isnan(nLogLs(iPtpnt,iRep))) ...
                || (any(~isnan(pars(iPtpnt, :, iRep))))
            
            if (~overrideDupError) || strcmp(firstOrAll, 'first')
                error('There already seems to be data here.')
            elseif overrideDupError
                warning(['Duplicates found. Overwriting rep number.', ...
                    ' Possible bug.'])
                while (~isnan(nLogLs(iPtpnt,iRep))) ...
                    || (any(~isnan(pars(iPtpnt, :, iRep))))
                    iRep = iRep + 1;
                    [pars, nLogLs] = expandArrays(pars, nLogLs, iRep, ...
                        Nptpnts, Npars);
                end
            else
                error('Bug')
            end
        end
        nLogLs(iPtpnt,iRep) = f.nLogL;
        pars(iPtpnt,:,iRep) = f.pars;
    end
    
    save([shortName,'.mat'],'pars','nLogLs')
end

if strcmp(firstOrAll, 'first')
    assert(~any(isnan(pars(:))))
    assert(~any(isnan(nLogLs(:))))
end

end

function [pars, nLogLs] = expandArrays(pars, nLogLs, iRep, Nptpnts, Npars)
% Do we need to expand the pars and nLogLs arrays?

if (iRep > size(pars, 3))
    extraRows = iRep - size(pars, 3);
    pars = cat(3, pars, nan(Nptpnts, Npars, extraRows));
    nLogLs = cat(2, nLogLs, nan(Nptpnts, extraRows));
end

end

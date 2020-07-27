function accumulate_pars(dataDir)

[~, Nptpnts] = getData(dataDir);
Config = load('Config.mat');
Nreps = Config.Nreps;

for itype = 1 : length(Config.ModelList)
    fname = ['pars/pars_' Config.ModelList{itype}];
    files = dir([fname,'_*']);
    
    % Load the first file just to find out the number of parameters
    f = load(fullfile(files(1).folder,files(1).name));
    Npars = length(f.pars);
    pars = nan(Nptpnts,Npars,Nreps);
    nLogLs = nan(Nptpnts,Nreps);
    
    for iFile = 1:length(files)
        f = load(fullfile(files(iFile).folder,files(iFile).name));
        fparts = split(files(iFile).name,{'_','.'});
        iSubj = str2double(fparts{end-2});
        iRep = str2double(fparts{end-1});
        if isfield(f, 'likelihood')
            nLogLs(iSubj,iRep) = f.likelihood;
            % TODO Remove this option. 
        else
            nLogLs(iSubj,iRep) = f.nLogL;
        end
        pars(iSubj,:,iRep) = f.pars;
    end
    
    save([fname,'.mat'],'pars','nLogLs')
end
function accumulate_pars

Nsubjs = 9;
Nreps = 20;
Npars = 4;

for itype = 1:2
    switch itype
        case 1
            fname = 'pars/pars_ibs_Bayes';
        case 2
            fname = 'pars/pars_ibs_Freq';
    end
    files = dir([fname,'_*']);
    pars = nan(Nsubjs,Npars,Nreps);
    likelihoods = nan(Nsubjs,Nreps);
    for iFile = 1:length(files)
        f = load(fullfile(files(iFile).folder,files(iFile).name));
        fparts = split(files(iFile).name,{'_','.'});
        iSubj = str2double(fparts{end-2});
        iRep = str2double(fparts{end-1});
        likelihoods(iSubj,iRep) = f.likelihood;
        pars(iSubj,:,iRep) = f.pars;
    end
    save([fname,'.mat'],'pars','likelihoods')
end
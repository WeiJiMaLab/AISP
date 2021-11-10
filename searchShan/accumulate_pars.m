function accumulate_pars

Nsubjs = 9;

for itype = 1:4
    switch itype
        case 1
            fname = 'pars/pars_ibs_Bayes';
            Npars = 4;
        case 2
            fname = 'pars/pars_ibs_Freq';
            Npars = 4;
        case 3
            fname = 'pars/pars_ibs_Sample';
            Npars = 5;
        case 4
            fname = 'pars/pars_ibs_cssample';
            Npars = 5;
    end
    files = dir([fname,'_*']);
    Nreps = 1;
    for iFile = 1:length(files)
        fparts = split(files(iFile).name,{'_','.'});
        iRep = str2double(fparts{end-1});
        if iRep > Nreps
            Nreps = iRep;
        end
    end
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
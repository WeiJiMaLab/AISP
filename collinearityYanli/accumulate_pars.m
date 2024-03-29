function accumulate_pars

Nsubjs = 8;

for itype = 1:5
    Nreps = 1;
    switch itype
        case 1
            fname = 'pars/pars_Bayes';
            Npars = 7;
        case 2
            fname = 'pars/pars_Freq';
            Npars = 7;
        case 3
            fname = 'pars/pars_Freq2';
            Npars = 7;
        case 4
            fname = 'pars/pars_sample';
            Npars = 8;
        case 5
            fname = 'pars/pars_cssample';
            Npars = 8;
    end
    files = dir([fname,'_*']);
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
        if f.likelihood == 0
            keyboard
        end
    end
    save([fname,'.mat'],'pars','likelihoods')
end
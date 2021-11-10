function accumulate_pars

Nsubjs = 8;
Npars = 7;

Nreps = 1;
for itype = 1:3
    switch itype
        case 1
            fname = 'pars/pars_Bayes';
        case 2
            fname = 'pars/pars_Freq';
        case 3
            fname = 'pars/pars_Freq2';
        case 3
            fname = 'pars/pars_sample';
        case 3
            fname = 'pars/pars_cssample';
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
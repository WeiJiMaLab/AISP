function accumulate_pars

Nsubjs = 11;
Nreps = 10;
Npars = 9;

for itype = 1:10
    switch itype
        case 1
            fname = 'pars/parsBayes';
        case 2
            fname = 'pars/parsFreq';
        case 3
            fname = 'pars/parsFreq2';
        case 4
            fname = 'pars/parsFreq3';
        case 5
            fname = 'pars/parsVar';
        case 6
            fname = 'pars/pars_ibs_Bayes';
        case 7
            fname = 'pars/pars_ibs_Freq';
        case 8
            fname = 'pars/pars_ibs_Freq2';
        case 9
            fname = 'pars/pars_ibs_Freq3';
        case 10
            fname = 'pars/pars_ibs_Var';
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
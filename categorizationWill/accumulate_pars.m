function accumulate_pars

Nsubjs = 11;
Nreps = 10;
Npars = 9;

for itype = 1:4
    switch itype
        case 1
            fname = 'pars/parsBayes';
        case 2
            fname = 'pars/parsFreq';
        case 3
            fname = 'pars/parsFreq2';
        case 4
            fname = 'pars/parsFreq3';
    end
    files = dir([fname,'_*']);
    pars = zeros(Nsubjs,Npars,Nreps);
    likelihoods = zeros(Nsubjs,Nreps);
    for iFile = 1:length(files)
        f = load(fullfile(files(iFile).folder,files(iFile).name));
        fparts = split(files(iFile).name,{'_','.'});
        iSubj = str2num(fparts{2});
        iRep = str2num(fparts{3});
        likelihoods(iSubj,iRep) = f.likelihood;
        pars(iSubj,:,iRep) = f.pars;
    end
    save([fname,'.mat'],'pars','likelihoods')
end
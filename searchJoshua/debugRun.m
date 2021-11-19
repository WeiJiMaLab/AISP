function debugRun(modelNum, type, runChecks)

addReqPaths()
Config = produceConfig();

Dirs = produceDirs(6);
[DSet, Nptpnts] = getData(Dirs.dataDir);


if strcmp(type, 'min')
    DatSubj = DSet.P(1).Data;
    designMat = struct2DesignMat(DatSubj, 'to matrix', true);
    pars = [0, 0, 0, 0, 0, 2, 0.1, 100];
    
    tic
    for i = 1 : 100
        aisp_simResponseWrapper(designMat, pars, ...
            Config.ModelList{modelNum}, runChecks);
    end
    toc

elseif strcmp(type, 'full')
    iRep = 1;
    iPtpnt = 1;
    idx = 1;
    fit_cluster_ibs(iRep, iPtpnt, Config.ModelList{modelNum}, DSet, ...
        idx, true, runChecks)

end



function debugRun(modelNum, type, runChecks)

addReqPaths()
Config = produceConfig();

Dirs = produceDirs(6);
[DSet, Nptpnts] = getData(Dirs.dataDir);


if strcmp(type, 'min')
    tic
    for i = 1 : 1000
        pars = [0, 0, 0, 0, 0, 2, 0.1, 10];
        ParamStruct = paramVec2Struct(pars, Config.ModelList{modelNum}, ...
            'to struct');
        aisp_simResponse(Config.ModelList{modelNum}, ParamStruct, ...
            DSet.P(1).Data, runChecks);
    end
    toc

elseif strcmp(type, 'full')
    iRep = 1;
    iPtpnt = 1;
    idx = 1;
    fit_cluster_ibs(iRep, iPtpnt, Config.ModelList{modelNum}, DSet, ...
        idx, true, runChecks)

end



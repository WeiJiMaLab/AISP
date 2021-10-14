classdef utils
    
    methods(Static)
        function [stimulus, response, performance] = readdata(subjid, half)
            % input argument half specifies which half of the data should
            % be returned, 0 means returning the whole data set, 1 means
            % returning the half with odd indices, 2 means returning the half
            % of data with even indices
            
            if ~exist('half','var')
                half = 0;
            end
            if ~ismember(half, [0,1,2])
                error('Input of half is not valid, please input 0,1,2.')
            end
            % get list of filenames
            if ismember(subjid, {'MBC','MG','RC','WYZ','XLM','YC','YL','YMH','YZ'})
                datapath = 'data/Experiment 1/';
            else
                datapath = 'data/Experiment 2/';
            end
            files = dir([datapath subjid '*.mat']);
            
            % read every file
            stimulus = [];
            response = [];
            performance = [];
            for ii=1:length(files)
                load([datapath files(ii).name]);
                stimulus = [stimulus; data(:,1:2)];
                response = [response; data(:,3)];
                performance = [performance; data(:,4)];
            end
            idx = length(stimulus);
            if half == 1
                stimulus = stimulus(1:2:idx-1,:);
                response = response(1:2:idx-1,:);
                performance = performance(1:2:idx-1,:);
            elseif half == 2
                stimulus = stimulus(2:2:idx,:);
                response = response(2:2:idx,:);
                performance = performance(2:2:idx,:);
            end
        
        end
        
        function [stimulus, response, performance] = readdata2(subjid,session,noFiles)
            % read data from session 1
            if ~exist('noFiles','var')
                noFiles = Inf;
            end
            % get list of filenames
            files = dir(['output2/session' num2str(session) '/' subjid '*.mat']);
            % read every file
            stimulus = [];
            response = [];
            performance = [];
            for ii=1:min(noFiles,length(files))
                load(['output2/session' num2str(session) '/' files(ii).name]);
                stimulus = [stimulus; data(:,1:2)];
                response = [response; data(:,3)];
                performance = [performance; data(:,4)];
            end
        end
        
        function [stimulus, response] = readfakedata(subjid)
            load(['fake_data/fix_s/' subjid '.mat']);

            stimulus = data(:,1:4);
            response = data(:,5);

        end
        
        function create_joblist(varargin)

            fid = fopen('joblist.txt','w');
            npars = length(varargin);
            % change cell to num    
            switch npars
                case 1
                    dim1 = varargin(1);
                    dim1 = [dim1{:}];
                    for ii = dim1;
                        fprintf(fid,'null,%d\n',ii);
                    end
                case 2
                    dim1 = varargin(1);
                    dim2 = varargin(2);
                    dim1 = [dim1{:}];
                    dim2 = [dim2{:}];
                    for ii = dim1
                        for jj = dim2
                            fprintf(fid,'null,%d,%d\n',ii,jj);
                        end
                    end
                case 3
                    dim1 = varargin(1);
                    dim2 = varargin(2);
                    dim3 = varargin(3);
                    dim1 = [dim1{:}];
                    dim2 = [dim2{:}];
                    dim3 = [dim3{:}];
                    for ii = dim1
                        for jj = dim2
                            for kk = dim3
                                fprintf(fid,'null,%d,%d,%d\n',ii,jj,kk);
                            end
                        end
                    end
                case 4
                    dim1 = varargin(1);
                    dim2 = varargin(2);
                    dim3 = varargin(3);
                    dim4 = varargin(4);
                    dim1 = [dim1{:}];
                    dim2 = [dim2{:}];
                    dim3 = [dim3{:}];
                    dim4 = [dim4{:}];
                    for ii = dim1
                        for jj = dim2
                            for kk = dim3
                                for mm = dim4
                                    fprintf(fid,'null,%d,%d,%d,%d\n',ii,jj,kk,mm);
                                end
                            end
                        end
                    end
                case 5
                    dim1 = varargin(1);
                    dim2 = varargin(2);
                    dim3 = varargin(3);
                    dim4 = varargin(4);
                    dim5 = varargin(5);
                    dim1 = [dim1{:}];
                    dim2 = [dim2{:}];
                    dim3 = [dim3{:}];
                    dim4 = [dim4{:}];
                    dim5 = [dim5{:}];
                    for ii = dim1
                        for jj = dim2
                            for kk = dim3
                                for mm = dim4
                                    for nn = dim5
                                        fprintf(fid,'null,%d,%d,%d,%d,%d\n',ii,jj,kk,mm,nn);
                                    end
                                end
                            end
                        end
                    end
            end
            fclose(fid);
        end
        
        function check_fit_pars(subj_idx,model_idx)
            % print out the parameters that touch the boundaries
            common_dir = 'real_data_results/fix_s/';
            common_dir2 = 'real_data_results/fix_s_0.001_0.3_31/';
            for ii = subj_idx
                subjid = compute.subjids{ii};
                for jj = model_idx
                    dirname = compute.model_names{jj};
                    % get parameters to generate the table for that subject and model
                    load([common_dir dirname '/pars/' subjid '_pars.mat'])
                    % get fit pars
                    load([common_dir dirname '/fit_pars/' subjid '.mat'])
                    if CPG_lambda_hat==max(lambda) || CPG_lambda_hat==lambda(1)
                        fprintf([subjid '_' dirname '_lambda:%6.3f\n'], CPG_lambda_hat)
                    end
                    
                end
            end
            
        end
            
        function show_fit_pars(subj_idx,model_idx,flex_s)
            if ~exist('flex_s','var')
                flex_s=0;
            end
            if flex_s
                common_dir = compute.dirs3;
            else
                common_dir = compute.dirs;
            end
            for jj = model_idx
                dirname = compute.model_names{jj};
                lambdaMat = zeros(1,length(subj_idx));
                guessMat = zeros(1,length(subj_idx));
                cnt = 0;
                for ii = subj_idx
                    subjid = compute.subjids{ii};
                    load([common_dir dirname '/fit_pars/' subjid '.mat']);
                    cnt = cnt + 1;
                    lambdaMat(cnt) = CPG_lambda_hat;
                    guessMat(cnt) = CPG_guess_hat;
                    if flex_s
                        fprintf([subjid ':%6.3f,%6.3f,%6.3f\n'], CPG_lambda_hat, CPG_std_s_hat, CPG_guess_hat);
                    else
                        fprintf([subjid ':%6.3f,%6.3f\n'], CPG_lambda_hat, CPG_guess_hat);
                    end
                end
                lambdaMean = mean(lambdaMat);
                lambdaStd = std(lambdaMat);
                guessMean = mean(guessMat);
                guessStd = std(guessMat);
                if flex_s
                    stdsMean = mean(CPG_std_s_hat);
                    stdsStd = std(CPG_std_s_hat);
                    fprintf([dirname ' lambda :%6.3f +- %6.3f\n'],lambdaMean, lambdaStd);
                    fprintf([dirname ' guess :%6.3f +-  %6.3f\n'],guessMean, guessStd);
                    fprintf([dirname ' sigmas :%6.3f +-  %6.3f\n'],stdsMean, stdsStd);
                else
                    fprintf([dirname ' lambda :%6.3f +- %6.3f\n'],lambdaMean, lambdaStd);
                    fprintf([dirname ' guess :%6.3f +-  %6.3f\n'],guessMean, guessStd);
                end
            end
        end
        
        function show_fit_pars_all
            parsMat = zeros(5,length(compute.subjids));
            for ii = 1:length(compute.subjids)
                subjid = compute.subjids{ii};
                dirname = 'flex2';
                load(['real_data_results/fix_s/' dirname '/fit_pars/' subjid '.mat']);
                parsMat(1,ii) = lambda_hat;
                parsMat(2,ii) = l1_hat;
                parsMat(3,ii) = l2_hat;
                parsMat(4,ii) = l3_hat;
                parsMat(5,ii) = guess_hat;
            end
            pars_mean = mean(parsMat,2);
            pars_std = std(parsMat,[],2)/sqrt(length(compute.subjids));
            
            pars_mean
            pars_std
            parsMat
        end
        
        function show_fake_pars(subj_idx, run_idx)
            % show parameters used to generate fake data
            for ii = subj_idx
                for jj = run_idx
                    subjid = [compute.model_names{ii} '_' num2str(jj)];
                    load([compute.dirs_fakedata1 subjid]);
                    fprintf([subjid ' lambda :%6.3f, %6.3f\n'],lambda, guess_rate);
                end
            end
        end
        
        function bins = generate_bins(nBins)
            nBins = nBins-1;
            if mod(nBins,2)==1
                error('Number of bins should be an odd number');
            end
            bins = norminv(linspace(0,1,nBins+2), 0,compute.s_std);
            bins = bins(2:end-1);
            bins_plus = bins(length(bins)/2+1:end);
            bins_minus = wrev(bins(1:length(bins)/2));
            
            bins_new = zeros(1, nBins+1);
            idx_mid = length(bins)/2+1;
            for ii = 1:length(bins)/2
                bins_new(idx_mid+ii) = bins_plus(ii)*2 - bins_new(idx_mid+ii-1);
                bins_new(idx_mid-ii) = bins_minus(ii)*2 - bins_new(idx_mid-ii+1);
            end
            bins = bins_new;
           
%             save('bins_new', 'bins');
        end
        
        function [cnt, cnt_r, p_right] = compute_cnts(stimuli, response, bins)
            % count the number of trials and right-reporting trials in each bin
            % 2015-11-23 SS
            target_stimuli = stimuli(:,1);
            dist_stimuli = stimuli(:,2);
            target_stimuli_right = target_stimuli(response==1);
            dist_stimuli_right = dist_stimuli(response==1);
            idx1 = interp1(bins,1:length(bins),target_stimuli, 'nearest','extrap');
            idx2 = interp1(bins,1:length(bins),dist_stimuli, 'nearest','extrap');
            idx1_right = interp1(bins,1:length(bins),target_stimuli_right, 'nearest','extrap');
            idx2_right = interp1(bins,1:length(bins),dist_stimuli_right, 'nearest','extrap');
            cnt_r = zeros(length(bins), length(bins));
            cnt = zeros(length(bins), length(bins));
            for ii = 1:length(bins)
                for jj = 1:length(bins)
                    cnt(ii,jj) = length(response(idx1==ii & idx2==jj));
                    temp = target_stimuli(idx1_right==ii & idx2_right==jj);
                    cnt_r(ii,jj) = length(temp);
                end
            end
            p_right = cnt_r./cnt;
        end
        
        function check_decision_rule(x, sigma, nSteps)
            if size(x) ~= [1,4]
                error('Input orietations are not valid');
            end
            if ~exist('sigma', 'var')
                sigma = 5;
            end
            if ~exist('nSteps', 'var')
                nSteps = 1000;
            end
            sT_plus = linspace(0,90,nSteps);
            sT_minus = linspace(-90,0,nSteps);
            sD = linspace(-90,90,nSteps);
            sigma_s = 9.0593;
            setsize = size(x,2);
            denom_L_num = zeros(nSteps, setsize);
            R_num = zeros(nSteps, setsize);
            nom_L_num = zeros(nSteps, setsize);
            denom_anl = zeros(1, setsize);
            nom_anl = zeros(1, setsize);
            
            idx = 1:size(x,2);
            for ii = 1:size(x,2)
                denom_L_num(:,ii) = exp(-(x(ii) - sT_minus).^2/2/sigma^2).*exp(-sT_minus.^2/2/sigma_s^2);
                nom_L_num(:,ii) = exp(-(x(ii) - sT_plus).^2/2/sigma^2).*exp(-sT_plus.^2/2/sigma_s^2);
                sD_temp = repmat(sD',1,3);
                R_num(:,ii) = prod(exp(-bsxfun(@minus,x(idx~=ii),sD_temp).^2/2/sigma^2),2)'.*exp(-sD.^2/2/sigma_s^2);
                weight = exp(-x(ii)^2/2/(sigma^2+sigma_s^2))...
                    *exp(-(setsize-1)/2/sigma^2*(sum(x(idx~=ii).^2/(setsize-1)) - (1/(setsize-1)*sum(x(idx~=ii)))^2))...
                    *exp(-(1/(setsize-1)*sum(x(idx~=ii)))^2/2/(sigma^2/(setsize-1)+sigma_s^2));
                denom_anl(ii) = (.5 - .5*erf(x(ii)/sigma^2/sqrt(2*(1/sigma^2+1/sigma_s^2))))*weight;
                nom_anl(ii) = (.5 + .5*erf(x(ii)/sigma^2/sqrt(2*(1/sigma^2+1/sigma_s^2))))*weight;
            end
            ratio_num = sum(sum(nom_L_num).*sum(R_num))/sum(sum(denom_L_num).*sum(R_num));
            ratio_anl = sum(nom_anl)/sum(denom_anl);
            ratio_num
            ratio_anl
            
        end
        
        function generate_stimuli(nTrials)
            sT = normrnd(0,compute.s_std, 1,nTrials);
            sD = normrnd(0,compute.s_std, 1,nTrials);
            save('stimuli_small', 'sT', 'sD');
        end
        
        function xMat = generate_measurements(sigma)
            load('stimuli_small.mat');
            stimuli = [sT; repmat(sD,compute.setsize-1,1)];
            xMat = zeros(size(stimuli,1), size(stimuli,2),length(sigma));
            for ii = 1:length(sigma);
                xMat(:,:,ii) = normrnd(stimuli, sigma(ii));
            end
            save('measurements.mat', 'xMat', 'sigma');
        end
        
        function check_learning(subj_idx)
            nSessions = 3;
            perfMat = zeros(length(subj_idx),nSessions);
            for ii = subj_idx
                subjid = compute.subjids{ii};
                for jj = 1:nSessions
                    fname = ['output2/session' num2str(jj) '/' subjid '_' num2str(jj)];
                    load(fname);
                    performance = data(:,4);
                    perfMat(ii,jj) = mean(performance);
                end
            end
            fig = Figure(101,'size',[50,40]);
            if length(subj_idx)==1
                plot(1:nSessions, perfMat);
                title([subjid ' performance']);
                
            else
                colorvec = get(gca,'ColorOrder');
                colorvec(8,:)= [0.3,0.5,0.2];
                colorvec(9,:)= [0.9,0,0.4];
                colorvec2 = min(colorvec+0.4,1);
                meanPerf = mean(perfMat);
                
                for ii = 1:nSessions
%                     bar(ii,meanPerf(ii),'FaceColor',colorvec2(ii,:),'EdgeColor','None'); hold on
                    bar(ii,meanPerf(ii),'FaceColor',[0.5,0.5,0.5]); hold on
                end
                for ii = subj_idx
                    h = plot(1:nSessions,perfMat(ii,:),'Color',colorvec2(ii,:));
                    h2 = scatter(1:nSessions,perfMat(ii,:),'CData', colorvec(ii,:),'Marker','.');
                    set(h2,'SizeData',150);
                end
                set(gca,'XTick',[1,2,3]);  
                xlabel('Session number');
                ylabel('Performance');
                set(gca,'Clipping','off');
                ylim([0,1]);
            end
            [h,p]=ttest(perfMat(:,1)-perfMat(:,2))
            fig.cleanup; fig.save([plots.dirs_save3 'performance.eps']);
        end
        
        function test_RMSE(model_idx)
            for ii = model_idx
                load([compute.dirs 'opt/errMat.mat']);
                errMat_opt = errMat;
                dirname = compute.model_names{ii};
                load([compute.dirs dirname '/errMat.mat'],'errMat');
                [h,p] = ttest(errMat_opt - errMat)
            end
        end
        
        function check_LL_expression(x, sigma, sigma_s)
                        
            weight = zeros(size(x));
            N = length(x);
            idx = 1:N;
            
            % analytical
            term = zeros(1,N);
            for ii = 1:length(weight)
                xT = x(idx==ii); xD = x(idx~=ii);
                weight(ii) = xT^2/(sigma_s^2+sigma^2) + mean(xD)^2/(sigma_s^2+sigma^2/(N-1)) + var(xD,1)*(N-1)/sigma^2;
                term(ii) = 1+erf(xT/sigma^2/sqrt(2*(1/sigma^2+1/sigma_s^2)));
            end
            
            prefactor = 0.5*sqrt(2*pi/(1/sigma^2+1/sigma_s^2))*sqrt(2*pi/((N-1)/sigma^2+1/sigma_s^2));
            analytical = sum(term.*exp(-weight/2))*prefactor
            
            sTMat = linspace(0,100,10000); sDMat = linspace(-100,100,10000);
            % numerical
            term1 = zeros(1,N); term2 = zeros(1,N);
            for ii = 1:N
                xT = x(idx==ii); xD = x(idx~=ii);
                term1(ii) = sum(exp(-(xT-sTMat).^2/2/sigma^2).*exp(-sTMat.^2/2/sigma_s^2))*range(sTMat)/length(sTMat);
                temp = 1;
                for jj = 1:N-1
                    temp = temp.*exp(-(xD(jj)-sDMat).^2/2/sigma^2);
                end
                term2(ii) = sum(temp.*exp(-sDMat.^2/2/sigma_s^2)*range(sDMat)/length(sDMat));
                term(ii) = term1(ii)*term2(ii);
            end
            
            numerical = sum(term)
            
            ratio = numerical/analytical
        end
        
        function residual_analysis(subj_idx)
           
            for ii = 1:length(subj_idx)
                subjid = compute.subjids{subj_idx(ii)};
                load([compute.dirs 'opt/fit_results_2d2/' subjid '.mat']);
                errMat = prediction-p_right;
                fig = Figure(101,'size',[60,40]);
                imagesc(errMat); caxis([-0.1,0.1]); colorbar

                % compute p values of each data point (chi2 test)
                
                chi2stat = cntVec.*errMat.^2./prediction + cntVec.*errMat.^2./(1-prediction);
                p = 1-chi2cdf(chi2stat,1);
                
                % test whether p is an uniform distribution
                test_cdf = makedist('Uniform');
                [h,p2] = kstest(p(:),'cdf',test_cdf)
                
                title({[subjid ' residual'], ['p = ' num2str(p2,'%.02f')]})
                fig.cleanup
                dir = 'plots/';
                if ~exist(dir,'dir')
                    mkdir(dir);
                end
                box on
                fig.save([dir subjid '_residual.eps'])
            end
            

        end
        
        function deviation_chi2_test(subj_idx,model_idx)
            % chi2_test for deviations
            % 11-22-2015 SS
            for ii = 1:length(subj_idx)
                for jj = 1:length(model_idx)
                    dirname = compute.model_names{jj};
                    subjid = compute.subjids{ii};
                    load([compute.dirs dirname '/fit_results_2d2/' subjid]);
                    p_right(p_right==0) = 0.0001;
                    p_right(p_right==1) = 1-0.0001;
                    D_vec = cntVec.*(p_right.*log(p_right./prediction)+(1-p_right).*log((1-p_right)./(1-prediction)));
                    D = sum(D_vec(:));
                    p = 1-chi2cdf(D,80)
                end
            end
        end
        
        function [entropyvalue,entropyMat] = compute_G_entropy(cnt, cnt_r)
            % compute Grassberger entropy
            % 15-11-23 SS
            gamma = 0.577215;  
            entropyMat = zeros(size(cnt));
            for ii = 1:size(cnt,1)
                for jj = 1:size(cnt,2)
                    
                    if cnt_r(ii,jj) == 0 || cnt_r(ii,jj) == cnt(ii,jj)
                        entropyMat(ii,jj) = 0;
                        continue
                    end
                    G_N = -gamma - log(2) + sum(2./(2*(1:floor(cnt(ii,jj)/2))-1));
                    if cnt_r(ii,jj)==1
                        G_n1 = -gamma - log(2);
                    else
                        G_n1 = -gamma - log(2) + sum(2./(2*(1:floor(cnt_r(ii,jj)/2))-1));
                    end
                    if cnt_r(ii,jj)==cnt(ii,jj)-1
                        G_n2 = -gamma - log(2);
                    else
                        G_n2 = -gamma - log(2) + sum(2./(2*(1:floor((cnt(ii,jj)-cnt_r(ii,jj))/2))-1));
                    end
                    entropyMat(ii,jj) = G_N - (cnt_r(ii,jj)*G_n1+(cnt(ii,jj)-cnt_r(ii,jj))*G_n2)/cnt(ii,jj);
                end
            end
            entropyvalue = -sum(entropyMat(:).*cnt(:));
        end
        
        function [mle,mleMat] = compute_cross_entropy(cnt,cnt_r,prediction_hat)
            % compute cross entropy by computing the mle
            cnt_l = cnt - cnt_r;
            mleMat = log(prediction_hat).*cnt_r + log(1-prediction_hat).*cnt_l;
            mle = sum(mleMat(:));
        end
        
        function copy_files(model_idx)
            % copy tables from fake_data folders to main folder
            % 15-07-07 SS
            for ii = model_idx
                dirname = compute.model_names{ii};
                new_dir = [compute.dirs_fake1 dirname];
                if ~exist(new_dir,'dir')
                    mkdir(new_dir);
                end
                new_file = [new_dir '/tables.mat'];
                old_file = ['../fake_data_test2/fake_data_results/fix_s/' dirname '/tables.mat'];
                copyfile(old_file, new_file);
            end
        
        end
        
        function [flex_s, npars, nTrials, common_dir, subj_ids,run_idx] = parse_inputs(subj_idx, varargin)
            % from the variable argument inputs, parse what conditions the
            % analysis should be on
            % 15-07-08 SS
            
            % case of fix s and real data
            if isempty(varargin{1})
                flex_s = 0;
                npars = 2;
                nTrials = 2000;
                common_dir = compute.dirs;
                run_idx = 0;
                subj_ids = compute.subjids(subj_idx);
                return
            end
            
            % case of flex s and real data
            if ischar(varargin{1}) % first argument is 'flex' instead of runs of fake data
                flex_s = 1;
                npars = 3;
                nTrials = 2000;
                common_dir = compute.dirs3;
                run_idx = 0;
                subj_ids = compute.subjids(subj_idx);
                return
            end
            
            % case of flex s and fake data
            if length(varargin)==2
                flex_s = 1;
                npars = 3;
                nTrials = 3000;
                common_dir = compute.dirs_fake2;
                runs = varargin{1};
                runs = runs{:};
                run_idx = runs;
                subj_ids = utils.create_subjids(subj_idx,runs);
                return
            end
            
            % case of fix s and fake data
            flex_s = 0;
            npars = 2;
            nTrials = 3000;
            common_dir = compute.dirs_fake1;
            runs = varargin{1};
            runs = runs{:};
            run_idx = runs;
            subj_ids = utils.create_subjids(subj_idx,runs);

        end
        
        function subj_ids = create_subjids(subj_idx, runs)
            % create names of fake data sets
            % 15-07-08 SS
            subj_ids = cell(length(subj_idx),length(runs));
            for ii = 1:length(subj_idx)
                for jj = 1:length(runs)
                    subj_ids{ii,jj} = [compute.model_names{subj_idx(ii)} '_' num2str(runs(jj))];
                end
            end
            subj_ids = reshape(subj_ids,1,numel(subj_ids));
        end
        
        function [eviMat,models] = get_evi(subj_idx,model_idx,type,norm, varargin)
            % get the matrix for evidence, for certain subjects and model_idx
            % if the data are real data, returned mat has size length(subj_idx)xlength(model_idx)
            % if the data are fake data, returned mat has size length(subj_idx)xlength(run_idx)xlength(model_idx)
            % eviMat for real data is not normalized
            % eviMat for fake data is normalized by the model to generate fake data
            % output models are the names of models shown on axis labels
            % SS 2015-07-16
            [~,~,~,common_dir,~,run_idx] = utils.parse_inputs(subj_idx,varargin);
            cmp_type = plots.cmp{type};
            models = cell(1,length(model_idx));
            for ii = 1:length(model_idx)
                models{ii} = plots.model_names{model_idx(ii)};
            end
                
            if strcmp(common_dir, compute.dirs)
                % real data
                eviMat = zeros(length(subj_idx),length(model_idx));
                for ii = 1:length(subj_idx)
                    for jj = 1:length(model_idx)
                        subjid = compute.subjids{subj_idx(ii)};
                        dirname = compute.model_names{model_idx(jj)};
                        load([common_dir dirname '/evi/' subjid]);
                        if type==3
                            eviMat(ii,jj) = eval(cmp_type);
                        else
                            eviMat(ii,jj) = -eval(cmp_type);
                        end
                    end
                end
            else
                % fake data
                eviMat = zeros(length(subj_idx),length(run_idx),length(model_idx));
                for ii = 1:length(model_idx)
                    for jj = 1:length(run_idx)
                        subjid = [compute.model_names{model_idx(ii)} '_' num2str(run_idx(jj))];
                        for kk = 1:length(model_idx)
                            dirname = compute.model_names{model_idx(kk)};
                            load([common_dir dirname '/evi/' subjid]);
                            eviMat(ii,jj,kk) = eval(cmp_type);
                        end
                        if norm
                            eviMat(ii,jj,:) = eviMat(ii,jj,:) - eviMat(ii,jj,ii);
                        end
                    end
                end % end of for model_idx
            end
         
        end
        
        function [lambda, guess_rate] = get_fit_pars(subj_idx, model_idx)
            % now only works for fix s and real data, needs to be
            % generalized to flex s and fake data in the future
            % SS 2015-07-16
            
            lambda = zeros(length(subj_idx),length(model_idx));
            guess_rate = zeros(size(lambda));
            
            for ii = 1:length(subj_idx)
                subjid = compute.subjids{subj_idx(ii)};
                for jj = 1:length(model_idx)
                    dirname = compute.model_names{model_idx(jj)};
                    load([compute.dirs dirname '/fit_pars/' subjid '.mat'])
                    lambda(ii,jj) = CPG_lambda_hat;
                    guess_rate(ii,jj) = CPG_guess_hat;
                end
            end
        end  
    
        function [lambda, guess_rate] = generate_fake_pars(subj_idx)
            % generate parameters based on parameters of real subjects
            % parameters are real pars of different models weighted by
            % their likelihood
            % SS 2015-07-16
            
            % load parameters and evidence
            [lambdaMat, guess_rateMat] = utils.get_fit_pars(subj_idx,1:25);
            eviMat = utils.get_evi(subj_idx,1:25,3);
            eviMat = bsxfun(@minus, eviMat, max(eviMat,[],2));
            lambda = sum(lambdaMat.*exp(eviMat),2)./sum(exp(eviMat),2);
            lambda = lambda';
            guess_rate = sum(guess_rateMat.*exp(eviMat),2)./sum(exp(eviMat),2);
            guess_rate = guess_rate';
            
        end
       
    end
    
end


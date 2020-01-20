classdef plots
    properties(Constant)
        % model comparison
        cmp = {'aic', 'bic', 'bmc','aicc'};
        % shared directory
        dirs = 'real_data_results/';
        model_names = {'Opt','Sum','Max', 'Min', 'Var', 'MaxT2', 'MaxT12', 'MaxT13','MaxT23','MaxT123',...
            'SumErf', 'SumErfT1','SumErfT2','SumErfT3','SumErfT12', 'SumErfT13','SumErfT23',...
            'SumXT1','SumXT2','SumXT3','SumXT12','SumXT13','SumXT23','SumXT123','Sign','Random','Err','Sample'};
        dirs_save = '/users/shanshen/Dropbox/Shan_VR/Optimality/Results and figures/fake data/model similarity analysis/';
        dirs_save2 = '/users/shanshen/Dropbox/Shan_VR/Optimality/Results and figures/real data/2_sessions/';
        dirs_save3 = '/users/shanshen/Dropbox/Shan_VR/Optimality/Results and figures/real data/3_sessions/';
        dirs_save4 = '/Users/shanshen/Google Drive/Shan Optimality/Manuscripts/201507 Psychological Review Revision/separated figures_full version/raw figures/';
    end
    methods(Static)
        %% the following functions plot the results of model comparison
        function model_comparison_all_subjs(subj_idx, model_idx, type, tosave)
            % with statistics
            % load data
            nSubjs = length(subj_idx);
            nModels = length(model_idx);
 
            [eviMat,models] = utils.get_evi(subj_idx,model_idx,type);

            % plot results
            fig = Figure(101,'size',[150,80]);
            
            subplot(1,2,2);
            mat = eviMat';
            mat = bsxfun(@minus, mat, mat(1,:));
            if length(model_idx)==2
                h1 = bar(mat(end,:));
                set(h1, 'FaceColor',[0.3,0.3,1],'LineStyle','None');
            else
                barh(reshape(wrev(mat(:)),size(eviMat,2),size(eviMat,1)));
                set(gca,'Ytick',1:nModels);
                set(gca,'YTickLabel',wrev(models));
                xlim([-600,100]);
                title('Individual subjects');
                xlabel('Model log likelihood relative to optimal model')
            end
            subplot(1,2,1);
            mat = eviMat;
            mat = bsxfun(@minus, mat, mat(:,1));
           
            mean_mat = squeeze(mean(mat,1))
            
            if sum(model_idx==1)
                h = zeros(1,length(model_idx));
                p = zeros(1,length(model_idx));
                t = zeros(1,length(model_idx));
                for ii = 1:length(model_idx)
                    [h(ii),p(ii),~,stats] = ttest(eviMat(:,1)- eviMat(:,ii));
                    t(ii) = stats.tstat;
                end
            end

            err_mat = squeeze(std(mat,1)/sqrt(nSubjs))
            
            if length(model_idx)==2
                h = bar(mean_mat(2)); hold all
                set(h,'FaceColor',[0.7,0.7,0.7],'BarWidth',0.3);
                errorbar(mean_mat(2),err_mat(2),'k');
            else
                h = barh(wrev(mean_mat)); hold all;
                set(h,'FaceColor', [0.7,0.7,0.7]);
                g = herrorbar(wrev(mean_mat),1:length(mean_mat), wrev(err_mat),'k');
                set(gca,'Ytick',1:nModels);
                set(gca,'YTickLabel',wrev(models));
                xlim([-400,100]);
                xlabel('Model log likelihood relative to optimal model')
            end
            fig.cleanup;
            if tosave
                fig.save('ic.eps');
            end
        end
        function model_comparison_fake(model_idx,run_idx,type)
            %% matrix to store the data
            eviMat = zeros(length(model_idx),length(run_idx),length(model_idx));
            cmp_type = plots.cmp{type};
            models = cell(1,length(model_idx));
            %% load data
            cnt2 = 0;
            winMat = ones(length(model_idx),length(run_idx));
            for ii = 1:length(model_idx)
               
                for jj = 1:length(run_idx)
                    subjid = [compute.model_names{model_idx(ii)} '_' num2str(run_idx(jj))];
                    if ii==jj==1
                        models = cell(1,length(model_idx)); cnt = 1;
                    end
                    for kk = 1:length(model_idx)
                        dirname = compute.model_names{model_idx(kk)};
                        if ii==jj==1
                            models{cnt} = plots.model_names{model_idx(kk)};
                            cnt = cnt + 1;
                        end
                        load([compute.dirs_fake1 dirname '/evi/' subjid]);
                        eviMat(ii,jj,kk) = eval(cmp_type);
                    end
                    eviMat(ii,jj,:) = eviMat(ii,jj,:) - eviMat(ii,jj,ii);
                    if sum(eviMat(ii,jj,:)>0)
                        cnt2 = cnt2+1;
                        winMat(ii,jj)=0;
                    end
                end
            end % end of for model_idx
            %% plot results
            cnt2
            save(['win_' num2str(type)],'winMat')
            eviMatmean = squeeze(mean(eviMat,2));
            mean(eviMatmean(2:end,1))
            min(eviMatmean(2:end,1))
            max(eviMatmean(2:end,1))
            eviMatstd = squeeze(std(eviMat,[],2))/sqrt(length(run_idx));
            
            eviMatplot = max(eviMatmean, -100);
           
            fig = Figure(205,'size',[150,120]); imagesc(eviMatplot);colorbar
            set(gca, 'XTick',1:length(model_idx));
            set(gca, 'YTick',1:length(model_idx));
            set(gca, 'XTickLabel', models);
            set(gca, 'YTickLabel', models);
            eviMatmean
            eviMatstd
            fig.cleanup
            box on
            fig.save('confusion_matrix.eps')
            Mat = cell(length(model_idx),length(model_idx));
            eviMatmean(eviMatstd>=10) = round(eviMatmean(eviMatstd>=10));
            eviMatstd(eviMatstd>=10) = round(eviMatstd(eviMatstd>=10));
            
            for ii = 1:length(model_idx)
                for jj = 1:length(model_idx)
                    Mat(ii,jj) = cellstr([num2str(eviMatmean(ii,jj)) char(177) num2str(eviMatstd(ii,jj))]);
                end
            end
%             xlswrite([cmp_type '_fake_data.xls'],Mat);
        end % end of function model_comparison_all
        function model_comparison_all_subjs_all_models(type)
            cmp_type = plots.cmp{type};
            nModels = 24+19;
            models = cell(1,nModels);
            nSubjs = length(compute.subjids);
            eviMat = zeros(nSubjs,nModels);
            cnt = 0;
            for ii = 1:2*24
                if ii<=24
                    common_dir = compute.dirs; dirname = compute.model_names{ii};
                    cnt = cnt+1; models{cnt} = dirname;
                    disp(['ii=' num2str(ii) ' cnt=' num2str(cnt)]);
                elseif ~isempty(find(compute.invalid_idx==ii-24))
                    disp('This model is invalid for this computation!');
                    continue
                else
                    common_dir = compute.dirs3; dirname = compute.model_names{ii-24};
                    cnt = cnt+1; models{cnt} = [dirname '_flex_s'];
                    disp(['ii=' num2str(ii) ' cnt=' num2str(cnt)]);
                end
                for jj = 1:nSubjs
                    subjid = compute.subjids{jj};
                    load([common_dir dirname '/evi/' subjid]);
                    eviMat(jj,cnt) = eval(cmp_type);
                end
                
            end
            figure;
            set(gcf,'Position',get(gcf,'Position').*[.1 .1 1.5 1]);
            subplot(1,2,1);
            mat = eviMat';
            mat = bsxfun(@minus, mat, mat(1,:));
            barh(reshape(wrev(mat(:)),size(eviMat,2),size(eviMat,1)));
            set(gca,'Ytick',1:nModels);
            set(gca,'YTickLabel',wrev(models));
            title('Individual subjects');
            xlabel('Model log likelihood relative to optimal model')
            subplot(1,2,2);
            mat = eviMat;
            mat = bsxfun(@minus, mat, mat(:,1));
            mean_mat = squeeze(mean(mat,1));
            err_mat = squeeze(std(mat,1)/sqrt(nSubjs));
             barh(wrev(mean_mat)); hold all;
            herrorbar(wrev(mean_mat),1:length(mean_mat), wrev(err_mat));
            set(gca,'Ytick',1:nModels);
            set(gca,'YTickLabel',wrev(models));
            title('Subject average');
            xlabel('Model log likelihood in relative to optimal model')

        end
        function model_comparison_all_subjs_2(subj_idx, model_idx, type)
            %% load data
            nSubjs = length(subj_idx);
            nModels = length(model_idx);
            %% matrix to store the data
            eviMat = zeros(nSubjs, nModels);
            cmp_type = plots.cmp{type};
            models = cell(1,nModels);
%             model_names = {'opt','sum','max','min','corr','min_max','min_corr','max_corr','min_max_corr','err'};
            %% load data
            
            for ii = 1:nSubjs
                subjid = compute.subjids{subj_idx(ii)};
                load([compute.dirs 'random/evi/' subjid])
                evi_rand = eval(cmp_type);
                for jj = 1:nModels
                    
                    dirname = compute.model_names{model_idx(jj)};
                    models{jj} = dirname;
                    load([compute.dirs dirname '/evi/' subjid]);
                    eviMat(ii,jj) = eval(cmp_type)/evi_rand;
                end
                
            end
            
            %% plot results
            figure;
%             set(gcf,'Position',get(gcf,'Position').*[.1 .1 1.5 1]);
%             subplot(1,2,1);
%             mat = eviMat';
%             mat = bsxfun(@minus, mat, mat(1,:));
%             barh(reshape(wrev(mat(:)),size(eviMat,2),size(eviMat,1)));
%             set(gca,'Ytick',1:nModels);
%             set(gca,'YTickLabel',wrev(models));
%             title('Individual subjects');
%             xlabel('Model log likelihood in relative to optimal model')
%             subplot(1,2,2);
            mat = eviMat;
%             mat = bsxfun(@minus, mat, mat(:,1));
%            
            mean_mat = squeeze(mean(mat,1));
            
            err_mat = squeeze(std(mat,1)/sqrt(nSubjs));
            barh(wrev(mean_mat)); hold all;
            herrorbar(wrev(mean_mat),1:length(mean_mat), wrev(err_mat));
            set(gca,'Ytick',1:nModels);
            set(gca,'YTickLabel',wrev(models));
            title('Subject average');
            xlabel('Nomalized model log likelihood')

        end
        function model_comparison_pair(subj_idx,model_idx,type)
            % compare log model likelihoods of the same decision rule, with or without flex_s
            % subj_idx: 1:9, model_idx is a single number
            % SS 2014-03-17
            if numel(model_idx)~=1
                error('This plot function does not support the plotting results more than one model.');
            end
            if ~isempty(find(compute.invalid_idx==model_idx))
                error('This plot function does not support this model.')
            end
            cmp_type = plots.cmp{type};
            nSubjs = length(subj_idx);
            nModels = 2;
            models = cell(1,nModels);
            eviMat = zeros(nSubjs,nModels);
            eviOpt = zeros(1,nSubjs);
            for ii = 1:nModels
                dirname = compute.model_names{model_idx};
                if ii==1
                    common_dir = compute.dirs;
                    models{ii} = dirname;
                else
                    common_dir = compute.dirs3;
                    models{ii} = [dirname '_flex_s'];
                end
                
                for jj = 1:nSubjs
                    subjid = compute.subjids{jj};
                    load([common_dir dirname '/evi/' subjid]);
                    eviMat(jj,ii) = eval(cmp_type);
                    if ii==1
                    % load results of optimal model
                    load([compute.dirs 'opt/evi/' subjid]);
                    eviOpt(jj) = eval(cmp_type);
                    end
                end
                
            end
 
            figure;
            set(gcf,'Position',get(gcf,'Position').*[.1 .1 1.5 1]);
            subplot(1,2,1);
            mat = eviMat';
            mat = bsxfun(@minus, mat, eviOpt);
            barh(reshape(wrev(mat(:)),size(eviMat,2),size(eviMat,1)));
            set(gca,'Ytick',1:nModels);
            set(gca,'YTickLabel',wrev(models));
            title('Individual subjects');
            xlabel('Model log likelihood relative to optimal model')
            subplot(1,2,2);
            mat = eviMat;
            mat = bsxfun(@minus, mat, eviOpt');
            mean_mat = squeeze(mean(mat,1));
            err_mat = squeeze(std(mat,1)/sqrt(nSubjs));
            barh(wrev(mean_mat),0.3); hold all;
            herrorbar(wrev(mean_mat),1:length(mean_mat), wrev(err_mat));
            set(gca,'Ytick',1:nModels);
            set(gca,'YTickLabel',wrev(models));
            title('Subject average');
            xlabel('Model log likelihood in relative to optimal model')
        end
        function model_comparison_pair2(subj_idx,model_idx,type)
            % compare log model likelihoods of the same decision rule, with or without flex_s
            % subj_idx: 1:9, model_idx can be a vector
            % SS 2014-04-13

            model_idx = setdiff(model_idx, compute.invalid_idx);
            cmp_type = plots.cmp{type};
            nSubjs = length(subj_idx);
            nModels = length(model_idx);
            models = cell(1,nModels);
            eviMat = zeros(nSubjs,nModels);
            
            for ii = 1:nModels
                dirname = compute.model_names{model_idx(ii)};
                models{ii} = plots.model_names{model_idx(ii)};
                for jj = 1:nSubjs
                    subjid = compute.subjids{subj_idx(jj)};
                    load([compute.dirs dirname '/evi/' subjid]);
                    fix_evi = eval(cmp_type);
                    load([compute.dirs3 dirname '/evi/' subjid]);
                    flex_evi = eval(cmp_type);
                    eviMat(jj,ii) = flex_evi - fix_evi;
                end
                
            end
 
            figure;
            set(gcf,'Position',get(gcf,'Position').*[.1 .1 1.5 1]);
            subplot(1,2,1);
            mat = eviMat';
            
            barh(reshape(wrev(mat(:)),size(eviMat,2),size(eviMat,1)));
            set(gca,'Ytick',1:nModels);
            set(gca,'YTickLabel',wrev(models));
            title('Individual subjects');
            xlabel('Model log likelihood relative to fix model')
            subplot(1,2,2);
            mat = eviMat;
            
            mean_mat = squeeze(mean(mat,1));
            err_mat = squeeze(std(mat,1)/sqrt(nSubjs));
            barh(wrev(mean_mat),0.3); hold all;
            herrorbar(wrev(mean_mat),1:length(mean_mat), wrev(err_mat));
            set(gca,'Ytick',1:nModels);
            set(gca,'YTickLabel',wrev(models));
            title('Subject average');
            xlabel('Model log likelihood in relative to fix model')
        end
        %% the following functions plot the psychometric curves, data or both
        function psymetric_curve_1d_all(subj_idx,model_idx,flex_s,tosave)
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            if ~exist('tosave','var')
                tosave = 0;
            end
            
            if flex_s
                common_dir = compute.dirs3;
            else
                common_dir = compute.dirs;
            end
            nSubjs = length(subj_idx);
            for mm = model_idx
                dirname = compute.model_names{mm};
                load([common_dir dirname '/fit_results_1d/' compute.subjids{subj_idx(1)}]);
                %% load and plot data
                predictionMat = zeros(nSubjs,length(bins_target),length(bins_dist));
                p_rightMat = zeros(nSubjs, length(bins_target), length(bins_dist));
                cntMat = zeros(nSubjs, length(bins_target), length(bins_dist));
                
                for ii = 1:nSubjs
                    subjid = compute.subjids{subj_idx(ii)};
                    load([common_dir dirname '/fit_results_1d/' subjid]);
                    predictionMat(ii,:,:) = prediction;
                    p_rightMat(ii,:,:) = p_right;
                    cntMat(ii,:,:) = cntVec;
                end
                
                % compute RMSE
                err = (p_rightMat - predictionMat).^2;
                err_temp = err.*cntMat;
                total_cnt = squeeze(sum(sum(cntMat,2),3));
                err_temp = bsxfun(@rdivide, err_temp, total_cnt);
                err_temp = sqrt(sum(sum(err_temp,2),3));
                                
                RMSE = mean(err_temp)
    
                p_right = squeeze(mean(p_rightMat));
                prediction = squeeze(mean(predictionMat));
                err_data = squeeze(std(p_rightMat)/sqrt(nSubjs));
                err_pred     = squeeze(std(predictionMat,1)/sqrt(nSubjs));
                upper_pred   = prediction + err_pred;
                lower_pred   = prediction - err_pred;

                % fig = Figure(101,'size',[25,22]);               
                fig = Figure(104,'size',[35,28]); 
               
                colorvec = get(gca, 'ColorOrder');
                colorvec = min(colorvec+.65,1);
                for jj = 1:length(bins_dist)
                    patch([bins_target wrev(bins_target)], [squeeze(upper_pred(:,jj))' wrev(squeeze(lower_pred(:,jj))')], colorvec(jj,:), 'Linestyle','None'); hold on;
                end
                h = errorbar(repmat(bins_target,length(bins_dist),1)', p_right, err_data,'o');  hold on;
                set(h,'MarkerSize',4);
%                 legend('distractor<-5','-5<distractor<5', 'distractor >5','Location','Northwest');
                xlim([min(bins_target)-2 max(bins_target)+2]); 
%                 xlabel('stimulus/degs');
%                 ylabel('proportion of reporting "right"');
%                 set(gca,'XTickLabel','');
%                 set(gca,'YTickLabel','');
%                 box off; title([dirname '-RMSE=' num2str(RMSE)]);
                fig.cleanup;
                if tosave
                    if ~exist([plots.dirs_save2 '1d/'], 'dir')
                        mkdir([plots.dirs_save2 '1d/']);
                    end
                    fig.save([plots.dirs_save2 '1d/' plots.model_names{mm} '.eps']);
                end
            end
        end
        function psymetric_curve_1d(subj_idx,model_idx,flex_s,tosave)
            if ~exist('tosave', 'var')
                tosave = 0;
            end
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            if flex_s
                common_dir = compute.dirs3;
            else
                common_dir = compute.dirs;
            end
            
            for mm = model_idx
                for nn = subj_idx
                    dirname = compute.model_names{mm};
                    subjid = compute.subjids{nn};
                    
                    %% load data
                    load([common_dir dirname '/fit_results_1d/' subjid]);
                                                       
                    % compute RMSE
                    err = (p_right - prediction).^2;
                    total_cnt = sum(cntVec(:));
                    err_data = err.*cntVec/total_cnt;
                    RMSE = sqrt(sum(err_data(:)));
                                      
                    a = figure; 
                    set(gcf, 'Position', get(gcf, 'Position').*[1,1,0.8,0.8]);

                    plot(bins_target,prediction); hold on;
                    
                    plot(bins_target, p_right,'o');  hold on;
%                     legend('distractor<-5','-5<distractor<5', 'distractor >5','Location','Northwest');
                    xlim([min(bins_target)-2 max(bins_target)+2]); ylim([0,1]);
                    set(gca, 'XTickLabel',[]);
                    set(gca, 'YTickLabel',[]);
                    if flex_s
                        title([subjid ', ' dirname ', flex sigma s']);
                    else
                        title([subjid ', ' dirname ', fix sigma s']);
                    end
                    xlabel('target orientation/degs');
                    ylabel('proportion of reporting "right"');
                    box off;RMSE
                    if tosave
                        saveas(a, ['/Users/sshan/My Documents/Dropbox/Shan_Variable_Precision/Optimality/raw plots for figures/individual/1d/psy_' dirname subjid '_1d'], 'emf');
                    end
                end
            end
        end      
        function psymetric_curve_1d_dist_all(subj_idx,model_idx,flex_s,tosave)
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            if ~exist('tosave','var')
                tosave = 0;
            end
            if flex_s
                common_dir = compute.dirs3;
            else
                common_dir = compute.dirs;
            end
            nSubjs = length(subj_idx);
            for mm = model_idx
                dirname = compute.model_names{mm};
                load([common_dir dirname '/fit_results_1d_dist/' compute.subjids{subj_idx(1)}]);
                %% load and plot data
                predictionMat = zeros(nSubjs,length(bins_target),length(bins_dist));
                p_rightMat = zeros(nSubjs, length(bins_target), length(bins_dist));
                cntMat = zeros(nSubjs, length(bins_target), length(bins_dist));
                
                for ii = 1:nSubjs
                    subjid = compute.subjids{subj_idx(ii)};
                    load([common_dir dirname '/fit_results_1d_dist/' subjid]);
                    predictionMat(ii,:,:) = prediction;
                    p_rightMat(ii,:,:) = p_right;
                    cntMat(ii,:,:) = cntVec;
                end
                
                % compute RMSE
                err = (p_rightMat - predictionMat).^2;
                err_temp = err.*cntMat;
                total_cnt = sum(sum(cntMat,1),3);
                err_temp = bsxfun(@rdivide, err_temp, total_cnt);
                err_temp = sqrt(sum(sum(err_temp,1),3));
                RMSE = mean(err_temp)
    
                p_right = squeeze(mean(p_rightMat));
                prediction = squeeze(mean(predictionMat));
                err_data = squeeze(std(p_rightMat)/sqrt(nSubjs));
                err_pred     = squeeze(std(predictionMat)/sqrt(nSubjs));
                upper_pred   = prediction + err_pred;
                lower_pred   = prediction - err_pred;

                                
                fig = Figure(106,'size',[35,28]);
                colorvec = get(gca, 'ColorOrder');
                colorvec = min(colorvec+.65,1);
                for jj = 1:length(bins_target)
                    patch([bins_dist wrev(bins_dist)],[squeeze(upper_pred(jj,:)) wrev(squeeze(lower_pred(jj,:)))],  colorvec(jj,:), 'Linestyle','None'); hold on;
                end
                h = errorbar(repmat(bins_dist,length(bins_target),1)',p_right', err_data','o');  hold on;
                set(h,'MarkerSize',4);
%                 legend('target<-5','-5<target<5', 'target>5','Location','Northwest');
                xlim([min(bins_dist)-2 max(bins_dist)+2]); 
%                 xlabel('stimulus/degs');
%                 ylabel('proportion of reporting "right"');
                fig.cleanup;
                if ~exist([plots.dirs_save2 '1d_dist/'], 'dir')
                    mkdir([plots.dirs_save2 '1d_dist/']);
                end
                if ~exist('tosave','var')
                    tosave = 0;
                end
                if tosave
                    fig.save([plots.dirs_save2 '1d_dist/' plots.model_names{mm} '.eps']);
                end
            end
        end
        function psymetric_curve_1d_dist(subj_idx,model_idx,flex_s,tosave)
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            if ~exist('tosave','var')
                tosave = 0;
            end
            if flex_s
                common_dir = compute.dirs3;
            else
                common_dir = compute.dirs;
            end
            for mm = model_idx
                for nn = subj_idx
                    dirname = compute.model_names{mm};
                    subjid = compute.subjids{nn};
                    load([common_dir dirname '/fit_results_1d_dist/' subjid]);
                    %% load and plot data
               
                   % compute RMSE
                    err = (p_right - prediction).^2;
                    err = err.*cntVec/sum(cntVec(:));
                    
                    RSE = sqrt(sum(err(:)))

                    a = figure; 
                    set(gcf, 'Position', get(gcf, 'Position').*[1,1,0.8,0.8]);

                    
                    plot(bins_dist, prediction); hold on;
                    
                    plot(bins_dist, p_right, 'o');  hold on;
    %                 legend('target<-5','-5<target<5', 'target>5','Location','Northwest');
                    xlim([min(bins_dist)-2 max(bins_dist)+2]); ylim([0,1]);
                    xlabel('distractor orientation/degs');
                    ylabel('proportion of reporting "right"');
                    set(gca,'XTickLabel','');
                    set(gca,'YTickLabel','');
                    if flex_s
                        title([subjid ', ' dirname ', flex sigma s']);
                    else
                        title([subjid ', ' dirname ', fix sigma s']);
                    end
                    box off;
                    if tosave
                        saveas(a, ['C:/Documents and Settings/sshan/My Documents/Dropbox/Shan_Variable_Precision/Optimality/raw plots for figures/individual/1d_dist/psy_' dirname '_' subjid '_1d_dist'], 'emf');
                    
                    end
                end
            end
        end
        function psymetric_curve_2d_all(model_idx)
            
            nSubjs= length(compute.subjids);
                
            for mm = model_idx
                dirname = compute.model_names{mm};
                load([compute.dirs dirname '/fit_results_2d/' compute.subjids{1}]);

                %% load and plot data
                predictionMat = zeros(nSubjs,length(bins),length(bins));
                p_rightMat = zeros(nSubjs, length(bins), length(bins));
                cntMat = zeros(nSubjs, length(bins), length(bins));
                for ii = 1:nSubjs
                    subjid = compute.subjids{ii};
                    load([compute.dirs dirname '/fit_results_2d/' subjid]);
                    predictionMat(ii,:,:) = prediction;
                    p_rightMat(ii,:,:) = p_right;
                    cntMat(ii,:,:) = cntVec;
                end
                err = (p_rightMat - predictionMat).^2;
                err_temp = err.*cntMat;
                total_cnt = squeeze(sum(sum(cntMat,2),3));
                err_temp = bsxfun(@rdivide, err_temp, total_cnt);
                err_temp = sqrt(sum(sum(err_temp,2),3));
                
                RMSE = mean(err_temp);

                p_right = squeeze(mean(p_rightMat));
                prediction = squeeze(mean(predictionMat));
               
                figure; 
                set(gcf, 'Position', get(gcf, 'Position').*[1,1,2,1]);

                a  = subplot(1,2,1);
                imagesc(bins, bins, p_right'); caxis([0,1]); colorbar
                set(gca, 'XTick', bins); set(gca, 'YTick', bins);
                set(gca, 'XTickLabel', bins); set(gca, 'YTickLabel', bins);
                xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
                title(a,'Psychometric function - data');
                b =  subplot (1,2,2);
                imagesc(bins, bins, prediction'); caxis([0,1]); colorbar
                set(gca, 'XTick', bins); set(gca, 'YTick', bins);
                set(gca, 'XTickLabel', bins); set(gca, 'YTickLabel', bins);
                xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
                title(b,[' Psychometric function model prediction - RMSE = ' num2str(RMSE)]);
               
            end
        end        
        function psymetric_curve_2d2_all(subj_idx,model_idx,flex_s,tosave)
            
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            
            if ~exist('tosave','var')
                tosave = 0;
            end
            
            if flex_s
                common_dir = compute.dirs3;
            else
                common_dir = compute.dirs;
            end
            
            
            nSubjs= length(subj_idx);
                
            for mm = model_idx
                dirname = compute.model_names{mm};
                load([common_dir dirname '/fit_results_2d2/' compute.subjids{1}]);
                %% load and plot data
                predictionMat = zeros(nSubjs,length(bins),length(bins));
                p_rightMat = zeros(nSubjs, length(bins), length(bins));
                cntMat = zeros(nSubjs, length(bins), length(bins));
                for ii = 1:nSubjs
                    subjid = compute.subjids{subj_idx(ii)};
                    load([common_dir dirname '/fit_results_2d2/' subjid]);
                    predictionMat(ii,:,:) = prediction;
                    p_rightMat(ii,:,:) = p_right;
                    cntMat(ii,:,:) = cntVec;
                end
                err = (p_rightMat - predictionMat).^2;
                err_temp = err.*cntMat;
                total_cnt = squeeze(sum(sum(cntMat,2),3));
                err_temp = bsxfun(@rdivide, err_temp, total_cnt);
                err_temp = sqrt(sum(sum(err_temp,2),3));
                errMat = err_temp
                save([common_dir dirname '/errMat.mat'],'errMat');
                RMSE = mean(err_temp)
                ste = std(err_temp)./sqrt(nSubjs)

                p_right = squeeze(mean(p_rightMat));
                prediction = squeeze(mean(predictionMat));
               
                fig=Figure(102,'size',[45,35]); 
               
%                 a  = subplot(1,2,1);
                bins_plot = -4:4;
                bins_names = cell(1,length(bins));
                for ii = 1:length(bins)
                    bins_names{ii} = sprintf('%10.1f',bins(ii));
                end
               
%                 imagesc(bins_plot, wrev(bins_plot), p_right'); caxis([0,1]); colorbar
%                 set(gca, 'XTick', bins_plot); set(gca, 'YTick',bins_plot);
%                 set(gca, 'XTickLabel', bins_names); set(gca, 'YTickLabel', wrev(bins_names));
%                 xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
%                 title(a,'Psychometric function - data');
%                 b =  subplot (1,2,2);
%                 imagesc(bins_plot, wrev(bins_plot), prediction'); caxis([0,1]); colorbar
%                 set(gca, 'XTick', bins_plot); set(gca, 'YTick', bins_plot);
%                 set(gca, 'XTickLabel', bins_names); set(gca, 'YTickLabel', wrev(bins_names));
%                 xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
%                 title(b,[' Psychometric function model prediction - RMSE = ' num2str(RMSE)]);
                
                imagesc(1:length(bins_plot), wrev(1:length(bins_plot)), prediction'); caxis([0,1]); colorbar
                set(gca, 'XTick', 1:length(bins_plot)); set(gca, 'YTick', 1:length(bins_plot));
                set(gca, 'XTickLabel', 1:length(bins_plot)); set(gca, 'YTickLabel', wrev(1:length(bins_plot)));
%                 xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
                title([ dirname ' - RMSE = ' num2str(RMSE)]);
                
                if ~exist([plots.dirs_save '2d/'], 'dir')
                    mkdir([plots.dirs_save '2d/']);
                end
                fig.cleanup
                if tosave
                    fig.save([plots.dirs_save2 '2d/' plots.model_names{mm}]);
                end
%                 figure; 
%                 set(gcf, 'Position', get(gcf, 'Position').*[1,1,2,0.8]);
% 
%                 a  = subplot(1,2,1);
%                 bins_plot = -4:4;
%                 bins_names = cell(1,length(bins));
%                 for ii = 1:length(bins)
%                     bins_names{ii} = sprintf('%10.1f',bins(ii));
%                 end
%                
%                 surf(bins_plot, wrev(bins_plot), p_right'); caxis([0,1]); colorbar
%                 set(gca, 'XTick', bins_plot); set(gca, 'YTick',bins_plot);
%                 set(gca, 'XTickLabel', bins_names); set(gca, 'YTickLabel', wrev(bins_names));
%                 xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
%                 title(a,'Psychometric function - data');
%                 b =  subplot (1,2,2);
%                 surf(bins_plot, wrev(bins_plot), prediction'); caxis([0,1]); colorbar
%                 set(gca, 'XTick', bins_plot); set(gca, 'YTick', bins_plot);
%                 set(gca, 'XTickLabel', bins_names); set(gca, 'YTickLabel', wrev(bins_names));
%                 xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
%                 title(b,[' Psychometric function model prediction - RMSE = ' num2str(RMSE)]);
               
            end
        end
        function psymetric_curve_2d2(model_idx, subj_idx)
                                      
            for mm = model_idx
                for nn = subj_idx
                    dirname = compute.model_names{mm};
                    subjid = compute.subjids{nn};
                    load([compute.dirs dirname '/fit_results_2d2/' subjid]);
                    err = (p_right - prediction).^2;
                    err = err.*cntVec/sum(cntVec(:));
                    RSE = sqrt(sum(err(:)))
                    a=figure; 
                    set(gcf, 'Position', get(gcf, 'Position').*[1,1,0.6,0.6]);

                    bins_plot = -4:4;
                    bins_names = cell(1,length(bins));
                    for ii = 1:length(bins)
                        bins_names{ii} = sprintf('%10.1f',bins(ii));
                    end


                    imagesc(1:length(bins_plot), wrev(1:length(bins_plot)), prediction'); caxis([0,1]);
                    set(gca, 'XTick', 1:length(bins_plot)); set(gca, 'YTick', 1:length(bins_plot));
                    set(gca, 'XTickLabel', 1:length(bins_plot)); set(gca, 'YTickLabel', wrev(1:length(bins_plot)));
    %                 xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
%                     title([ dirname ' - RMSE = ' num2str(RMSE)]);
%                     saveas(a, ['C:/Documents and Settings/sshan/My Documents/Dropbox/Shan_Variable_Precision/Optimality/raw plots for figures/individual/2d/psy_' dirname '_' subjid '_2d'], 'emf');

                end
            end
        end
        function performance_2d2_all(model_idx)
            
            nSubjs= length(compute.subjids);
                
            for mm = model_idx
                dirname = compute.model_names{mm};
                load([compute.dirs dirname '/fit_results_2d2/perf' compute.subjids{1}]);
                %% load and plot data
                predictionMat = zeros(nSubjs,length(bins),length(bins));
                p_rightMat = zeros(nSubjs, length(bins), length(bins));
                cntMat = zeros(nSubjs, length(bins), length(bins));
                for ii = 1:nSubjs
                    subjid = compute.subjids{ii};
                    load([compute.dirs dirname '/fit_results_2d2/perf' subjid]);
                    predictionMat(ii,:,:) = prediction;
                    p_rightMat(ii,:,:) = p_right;
                    cntMat(ii,:,:) = cntVec;
                end
                err = (p_rightMat - predictionMat).^2;
                err_temp = err.*cntMat;
                total_cnt = squeeze(sum(sum(cntMat,2),3));
                err_temp = bsxfun(@rdivide, err_temp, total_cnt);
                err_temp = sqrt(sum(sum(err_temp,2),3));
                RMSE = mean(err_temp);

                p_right = squeeze(mean(p_rightMat));
                prediction = squeeze(mean(predictionMat));
               
                a=figure; 
                set(gcf, 'Position', get(gcf, 'Position').*[1,1,0.6,0.6]);

%                 a  = subplot(1,2,1);
                bins_plot = -4:4;
                bins_names = cell(1,length(bins));
                for ii = 1:length(bins)
                    bins_names{ii} = sprintf('%10.1f',bins(ii));
                end
               
%                 imagesc(bins_plot, wrev(bins_plot), p_right'); caxis([0,1]); colorbar
%                 set(gca, 'XTick', bins_plot); set(gca, 'YTick',bins_plot);
%                 set(gca, 'XTickLabel', bins_names); set(gca, 'YTickLabel', wrev(bins_names));
%                 xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
%                 title(a,'Psychometric function - data');
%                 b =  subplot (1,2,2);
%                 imagesc(bins_plot, wrev(bins_plot), prediction'); caxis([0,1]); colorbar
%                 set(gca, 'XTick', bins_plot); set(gca, 'YTick', bins_plot);
%                 set(gca, 'XTickLabel', bins_names); set(gca, 'YTickLabel', wrev(bins_names));
%                 xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
%                 title(b,[' Psychometric function model prediction - RMSE = ' num2str(RMSE)]);
                
                imagesc(1:length(bins_plot), wrev(1:length(bins_plot)), prediction'); caxis([0,1]);
                set(gca, 'XTick', 1:length(bins_plot)); set(gca, 'YTick', 1:length(bins_plot));
                set(gca, 'XTickLabel', 1:length(bins_plot)); set(gca, 'YTickLabel', wrev(1:length(bins_plot)));
%                 xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
                title([ dirname ' - RMSE = ' num2str(RMSE)]);
                if ~exist([plots.dir_save 'performance/'],'dir')
                    mkdir([plots.dir_save 'performance/']);
                end
                saveas(a, [plots.dir_save 'perf_' dirname], 'emf');
                
%                 figure; 
%                 set(gcf, 'Position', get(gcf, 'Position').*[1,1,2,0.8]);
% 
%                 a  = subplot(1,2,1);
%                 bins_plot = -4:4;
%                 bins_names = cell(1,length(bins));
%                 for ii = 1:length(bins)
%                     bins_names{ii} = sprintf('%10.1f',bins(ii));
%                 end
%                
%                 surf(bins_plot, wrev(bins_plot), p_right'); caxis([0,1]); colorbar
%                 set(gca, 'XTick', bins_plot); set(gca, 'YTick',bins_plot);
%                 set(gca, 'XTickLabel', bins_names); set(gca, 'YTickLabel', wrev(bins_names));
%                 xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
%                 title(a,'Psychometric function - data');
%                 b =  subplot (1,2,2);
%                 surf(bins_plot, wrev(bins_plot), prediction'); caxis([0,1]); colorbar
%                 set(gca, 'XTick', bins_plot); set(gca, 'YTick', bins_plot);
%                 set(gca, 'XTickLabel', bins_names); set(gca, 'YTickLabel', wrev(bins_names));
%                 xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
%                 title(b,[' Psychometric function model prediction - RMSE = ' num2str(RMSE)]);
               
            end
        end
        function psymetric_curve_2d3_all(model_idx)
            
            nSubjs= length(compute.subjids);
                
            for mm = model_idx
                dirname = compute.model_names{mm};
                load([compute.dirs dirname '/fit_results_2d3/' compute.subjids{1}]);
                load('bins');
                %% load and plot data
                predictionMat = zeros(nSubjs,length(bins_fine),length(bins_fine));
                p_rightMat = zeros(nSubjs, length(bins), length(bins));
                
                for ii = 1:nSubjs
                    subjid = compute.subjids{ii};
                    load([compute.dirs dirname '/fit_results_2d3/' subjid]);
                    predictionMat(ii,:,:) = prediction;
                    p_rightMat(ii,:,:) = p_right;
                end
    

                p_right = squeeze(mean(p_rightMat));
                prediction = squeeze(mean(predictionMat));
               
                figure; 
                set(gcf, 'Position', get(gcf, 'Position').*[1,1,2,0.8]);

                a  = subplot(1,2,1);
                bins_plot = -4:4;
                bins_names = cell(1,length(bins));
                for ii = 1:length(bins)
                    bins_names{ii} = sprintf('%10.1f',bins(ii));
                end
               
                imagesc(bins_plot, wrev(bins_plot), p_right'); caxis([0,1]); colorbar
                set(gca, 'XTick', bins_plot); set(gca, 'YTick',bins_plot);
                set(gca, 'XTickLabel', bins_names); set(gca, 'YTickLabel', wrev(bins_names));
                xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
                title(a,'Psychometric function - data');
                b =  subplot (1,2,2);
                imagesc(bins_plot, wrev(bins_plot), prediction'); caxis([0,1]); colorbar
                set(gca, 'XTick', bins_plot); set(gca, 'YTick', bins_plot);
                set(gca, 'XTickLabel', bins_names); set(gca, 'YTickLabel', wrev(bins_names));
                xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
                title(b,' Psychometric function model prediction ');
               
            end
        end
        function psymetric_curve_1d_all_flex
            nSubjs = length(compute.subjids);
            
            dirname = 'flex2';
            load([compute.dirs dirname '/fit_results_1d/' compute.subjids{1}]);
            %% load and plot data
            predictionMat = zeros(nSubjs,length(bins_target),length(bins_dist));
            p_rightMat = zeros(nSubjs, length(bins_target), length(bins_dist));
            cntMat = zeros(nSubjs, length(bins_target), length(bins_dist));

            for ii = 1:nSubjs
                subjid = compute.subjids{ii};
                load([compute.dirs dirname '/fit_results_1d/' subjid]);
                predictionMat(ii,:,:) = prediction;
                p_rightMat(ii,:,:) = p_right;
                cntMat(ii,:,:) = cntVec;
            end

            % compute RMSE
            err = (p_rightMat - predictionMat).^2;
            err_temp = err.*cntMat;
            total_cnt = squeeze(sum(sum(cntMat,2),3));
            err_temp = bsxfun(@rdivide, err_temp, total_cnt);
            err_temp = sqrt(sum(sum(err_temp,2),3));
            RMSE = mean(err_temp);

            p_right = squeeze(mean(p_rightMat));
            prediction = squeeze(mean(predictionMat));
            err_data = squeeze(std(p_rightMat)/sqrt(nSubjs));
            err_pred     = squeeze(std(predictionMat,1)/sqrt(nSubjs));
            upper_pred   = prediction + err_pred;
            lower_pred   = prediction - err_pred;


            figure; 
            set(gcf, 'Position', get(gcf, 'Position').*[1,1,1,0.8]);

            colorvec = get(gca, 'ColorOrder');
            colorvec = min(colorvec+.65,1);
            for jj = 1:length(bins_dist)
                patch([bins_target wrev(bins_target)], [squeeze(upper_pred(:,jj))' wrev(squeeze(lower_pred(:,jj))')], colorvec(jj,:), 'Linestyle','None'); hold on;
            end
            errorbar(repmat(bins_target,length(bins_dist),1)', p_right, err_data,'o');  hold on;
            legend('distractor<-5','-5<distractor<5', 'distractor >5','Location','Northwest');
            xlim([min(bins_target) max(bins_target)]); 
            xlabel('stimulus/degs');
            ylabel('proportion of reporting "right"');
            box on; 
       
        end
        function psymetric_curve_1d_data(subj_idx)
            nSubjs = length(subj_idx);
            
            dirname = compute.model_names{1};
            load([compute.dirs dirname '/fit_results_1d/' compute.subjids{1}]);
            %% load and plot data
            p_rightMat = zeros(nSubjs, length(bins_target), length(bins_dist));
            
            for ii = 1:nSubjs
                subjid = compute.subjids{subj_idx(ii)};
                load([compute.dirs dirname '/fit_results_1d/' subjid]);
                p_rightMat(ii,:,:) = p_right;
            end

            % compute RMSE
            
            p_right = squeeze(mean(p_rightMat));
            p_right = p_right';
            err_data = squeeze(std(p_rightMat)/sqrt(nSubjs));
            
            fig = Figure(103,'size',[35,25]);

            plot(bins_target, p_right); hold on;
                        
            h = errorbar(repmat(bins_target,length(bins_dist),1)', p_right', err_data,'o');  hold on;
            set(h,'MarkerSize',4);
%             legend('distractor<-5','-5<distractor<5', 'distractor >5','Location','Northwest');
            xlim([min(bins_target)-2 max(bins_target)+2]); 
%             xlabel('stimulus/degs');
%             ylabel('proportion of reporting "right"');
            box off
            fig.cleanup; fig.save([plots.dirs_save2 '1d_data.eps']);
            
        end
        function psymetric_curve_1d_dist_data(subj_idx)
            nSubjs = length(subj_idx);
            
            dirname = compute.model_names{1};
            load([compute.dirs dirname '/fit_results_1d_dist/' compute.subjids{1}]);
            %% load and plot data
            p_rightMat = zeros(nSubjs, length(bins_target), length(bins_dist));
            
            for ii = 1:nSubjs
                subjid = compute.subjids{subj_idx(ii)};
                load([compute.dirs dirname '/fit_results_1d_dist/' subjid]);
                p_rightMat(ii,:,:) = p_right;
            end

            % compute RMSE
            
            p_right = squeeze(mean(p_rightMat));
            p_right = p_right';
            err_data = squeeze(std(p_rightMat)/sqrt(nSubjs));
            
            fig = Figure(105,'size',[35,28]);
           
            plot(bins_dist,p_right); hold on;
                        
            h = errorbar(repmat(bins_dist,length(bins_target),1)',p_right, err_data','o');  hold on;
            set(h,'MarkerSize',4);
%             legend('target<-5','-5<target<5', 'target>5','Location','Northwest');
            xlim([min(bins_dist)-2 max(bins_dist)+2]); ylim([0,1]);
            fig.cleanup; fig.save([plots.dirs_save2 '1d_dist_data.eps']);
%             xlabel('stimulus/degs');
%             ylabel('proportion of reporting "right"');
            
            
        end
        function psymetric_curve_2d2_data(subj_idx)
            nSubjs = length(subj_idx);
            subj_ids = compute.subjids(subj_idx);
            dirname = compute.model_names{1};
            load([compute.dirs dirname '/fit_results_2d2/' compute.subjids{1}]);
            %% load and plot data
            p_rightMat = zeros(nSubjs, length(bins), length(bins));
            for ii = 1:nSubjs
                subjid = subj_ids{ii};
                load([compute.dirs dirname '/fit_results_2d2/' subjid]);
                p_rightMat(ii,:,:) = p_right;
            end
            
            p_right = squeeze(mean(p_rightMat));
           
            fig = Figure(101,'size',[45,35]); 
            

            bins_plot = -4:4;
            bins_names = cell(1,length(bins));
            for ii = 1:length(bins)
                bins_names{ii} = sprintf('%10.1f',bins(ii));
            end
            imagesc(1:length(bins_plot), wrev(1:length(bins_plot)), p_right'); caxis([0,1]); colorbar
            set(gca, 'XTick', 1:length(bins_plot)); set(gca, 'YTick',1:length(bins_plot));
            set(gca, 'XTickLabel', 1:length(bins_plot)); set(gca, 'YTickLabel', wrev(1:length(bins_plot)));
%             xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
            title('Psychometric function - data');
            fig.cleanup; fig.save([plots.dirs_save2 '2d_data.eps'])
           
        end
        function psymetric_curve_2d_data(subj_idx)
            for mm = subj_idx
                subjid = compute.subjids{mm};
                dirname = compute.model_names{1};
                load([compute.dirs dirname '/fit_results_2d2/' subjid]);
                %% load and plot data

                a=figure; 
                set(gcf, 'Position', get(gcf, 'Position').*[1,1,.8,.8]);

                bins_plot = -4:4;
                bins_names = cell(1,length(bins));
                for ii = 1:length(bins)
                    bins_names{ii} = sprintf('%10.1f',bins(ii));
                end
                p_rightMat = flipud(p_right');
                imagesc(1:length(bins_plot), wrev(1:length(bins_plot)), p_rightMat); caxis([0,1]);
                save([subjid '_p_right.mat'], 'p_rightMat');
                set(gca, 'XTick', 1:length(bins_plot)); set(gca, 'YTick',1:length(bins_plot));
                set(gca, 'XTickLabel', 1:length(bins_plot)); set(gca, 'YTickLabel', wrev(1:length(bins_plot)));
    %             xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
%     saveas(a, ['C:/Documents and Settings/sshan/My Documents/Dropbox/Shan_Variable_Precision/Optimality/raw plots for figures/individual/data_' subjid '_2d'], 'emf');
                title('Psychometric function - data');
                
            end
        end
        function performance_curve_2d_data
            nSubjs= length(compute.subjids);
                
            dirname = compute.model_names{1};
            load([compute.dirs dirname '/fit_results_2d2/perf' compute.subjids{1}]);
            %% load and plot data
            p_rightMat = zeros(nSubjs, length(bins), length(bins));
            for ii = 1:nSubjs
                subjid = compute.subjids{ii};
                load([compute.dirs dirname '/fit_results_2d2/perf' subjid]);
                p_rightMat(ii,:,:) = p_right;
            end
            
            p_right = squeeze(mean(p_rightMat));
           
            figure; 
            set(gcf, 'Position', get(gcf, 'Position').*[1,1,.8,.8]);

            bins_plot = -4:4;
            bins_names = cell(1,length(bins));
            for ii = 1:length(bins)
                bins_names{ii} = sprintf('%10.1f',bins(ii));
            end
            imagesc(1:length(bins_plot), wrev(1:length(bins_plot)), p_right'); caxis([0,1]); 
            set(gca, 'XTick', 1:length(bins_plot)); set(gca, 'YTick',1:length(bins_plot));
            set(gca, 'XTickLabel', 1:length(bins_plot)); set(gca, 'YTickLabel', wrev(1:length(bins_plot)));
%             xlabel('Target orientation/degs'); ylabel('Distractor orientation/degs');
            title('Psychometric function - data');
           
        end
        function show_landscape_flex(model_idx,subj_idx)
            % show landscape of model log likelihood of parameters lambda
            % and lambda_s. Guessing rate is the estimated value.
            % SS 04-09-14
            for ii = model_idx
                if sum(compute.invalid_idx==ii);
                    disp('This model is invalid for this plotting function!')
                    continue
                end
                dirname = compute.model_names{ii};
                for jj = subj_idx
                    subjid = compute.subjids{jj};
                    % load LL table and parameters
                    load([compute.dirs3 dirname '/tables/' subjid '.mat']);
                    load([compute.dirs3 dirname '/pars/' subjid '_pars.mat']);
                    % load fit pars
                    load([compute.dirs3 dirname '/fit_pars/' subjid '.mat']);
                    idx_guess = guess==CPG_guess_hat;
                    temp = squeeze(log_CPG_prediction(:,:,idx_guess));
                    temp = temp - max(temp(:));
                    figure; set(gcf,'Position', get(gcf,'Position').*[1,1,1.2,0.5]);
                    subplot 121
                    imagesc(lambda, std_s_vec, temp'); caxis([-100,0]); colorbar
                    subplot 122
                    imagesc(lambda, std_s_vec, exp(temp')); caxis([exp(-100),1]); colorbar
                    title(['model likelihood ' subjid]);
                end
            end
            
        end
        function show_landscape_flex2(subj_idx)
            dirname = 'flex2';
            subjid = compute.subjids{subj_idx};
            load([compute.dirs dirname '/tables/' subjid]);
            LLH = log_CPG_prediction;
            LLH = LLH - max(LLH(:));
            [~,idx] = max(LLH(:));
            [idx1 idx2 idx3 idx4 idx5] = ind2sub(size(LLH), idx);
            l1 = flexible_model.lambda1_vec;
            l2 = flexible_model.lambda2_vec;
            l3 = flexible_model.lambda3_vec;
            figure;
            set(gcf,'Position',get(gcf,'Position').*[1,1,2,0.5]);
            subplot 131
            imagesc(l1,l2,squeeze(exp(LLH(idx1,:,:,idx4,idx5)))'); colorbar
            xlabel('l1'); ylabel('l2');
            subplot 132
            imagesc(l2,l3,squeeze(exp(LLH(idx1,idx2,:,:,idx5)))'); colorbar
            xlabel('l2'); ylabel('l3');
            subplot 133
            imagesc(l1,l3,squeeze(exp(LLH(idx1,:,idx3,:,idx5)))'); colorbar
            xlabel('l1'); ylabel('l3');
        end
        function show_complexity
            names = {'cal', 'sum', 'err', 'sort', 'exp'};
            Mat = zeros(10, length(names));
            Mat(:,1) = wrev([160, 0 , 0 , 0, 20, 16, 24, 40, 0, 120]/40);
            Mat(:,2) = wrev([4, 0, 0, 3, 0, 0, 0, 0, 1, 1]);
            Mat(:,3) = wrev([4, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
            Mat(:,4) = wrev([0, 1, 1, 1, 1, 1, 1, 1, 1, 0]);
            Mat(:,5) = wrev([12, 0, 0, 0, 2, 2, 2, 3, 0, 12]/3);
            nModels = 10;
            models = cell(1,nModels);
            for ii = 1:nModels
                models{ii} = compute.model_names{ii};
            end
            
            figure; set(gcf,'Position',get(gcf,'Position').*[.1 .1 1.5 1]);
            barh(Mat, 'stack');
            legend(names); 
            set(gca, 'YTickLabel', wrev(models));
            set(gca, 'XTickLabel', '');
            xlabel('Complexity'); ylabel('Models');
        end
        function stimuli_distribution
            x = linspace(-40,40,500);
            y = normpdf(x,0,compute.s_std);
            figure; set(gcf, 'Position' , get(gcf, 'Position').*[1,1,0.75,0.6]);
            plot(x,y, 'LineWidth',2); hold on;
            plot(x,y,'--r', 'LineWidth',3);
            set(gca,'FontSize',11, 'FontName', 'Arial');
            legend('Target', 'Distractor');
            xlabel('Target or distractor orientation', 'FontSize',11, 'FontName', 'Arial');
            ylabel('Probability density', 'FontSize',11, 'FontName', 'Arial');
        end
        function quantile
            x = linspace(-40,40,500);
            y = normpdf(x,0,compute.s_std);
            
            figure; set(gcf, 'Position' , get(gcf, 'Position').*[1,1,0.75,0.6]);
            plot(x,y, 'LineWidth',2); hold on;
            load('bins_plotting');
            y2 = normpdf(bins,0,compute.s_std);
            for ii = 1:length(bins)
                line([bins(ii),bins(ii)],[0,y2(ii)]); 
            end
            set(gca,'FontSize',11, 'FontName', 'Arial');
%             legend('Target', 'Distractor');
%             xlabel('Target or distractor orientation', 'FontSize',11, 'FontName', 'Arial');
%             ylabel('Probability density', 'FontSize',11, 'FontName', 'Arial');
        end
        function performance_curve(model_idx)
            sigma = 0;
            dirname_temp = compute.model_names{model_idx(1)};
            load([compute.dirs dirname_temp '/perf/performance.mat']);
            perfMat = zeros(length(model_idx), length(performance));
            for ii = 1:length(model_idx)
                dirname = compute.model_names{model_idx(ii)};
                load([compute.dirs dirname '/perf/performance.mat']);
                perfMat(ii,:) = performance;
            end
%             figure; plot(sigma,perfMat);
            figure_size = [40,28];
            fig1 = Figure(101,'size',figure_size);
            p1=plot(sigma, perfMat(1,:),'k');hold on;
            p2=plot(sigma, perfMat([2:5,25],:));
            set(p1,'LineWidth',0.5); set(p2,'LineWidth',0.5);
            ylim([0.5,1]); box off
            fig1.cleanup; fig1.save('performance_cat1');
%             legend('Opt', 'Sum', 'Max', 'Min', 'Var');
%             xlabel('sigma');
%             ylabel('similarity to opt model'); title('prediction similarity to opt model, category 1');
            fig2 = Figure(102,'size',figure_size);
            p1=plot(sigma, perfMat(1,:), 'k'); hold on;
            p2=plot(sigma, perfMat(6:10,:));
            set(p1,'LineWidth',0.5); set(p2,'LineWidth',0.5);
            ylim([0.5,1]); box off
            fig2.cleanup; fig2.save([plots.dirs_save 'performance_cat2']);
%             legend('Opt','MaxT2','MaxT12','MaxT13','MaxT23','MaxT123');
%             xlabel('sigma');
%             ylabel('similarity to opt model'); title('prediction similarity to opt model, category 2');
            fig3 = Figure(103,'size', figure_size);
            p1=plot(sigma, perfMat(1,:), 'k'); hold on;
            p2=plot(sigma, perfMat(11:16,:));
            p3=plot(sigma, perfMat(17,:), 'Color', [0.9,0.5,0.2]);
            set(p1,'LineWidth',0.5);set(p2,'LineWidth',0.5);set(p3,'LineWidth',0.5);
            ylim([0.5,1]); box off
            fig3.cleanup; fig3.save([plots.dirs_save 'performance_cat3']);
%             legend('Opt','SumErf','SumErfT1', 'SumErfT2','SumErfT3','SumErfT12','SumErfT13', 'SumErfT23');
%             xlabel('sigma');
%             ylabel('similarity to opt model'); title('prediction similarity to opt model, category 3');
            fig4 = Figure(104,'size',figure_size);
            p1=plot(sigma, perfMat(1,:), 'k'); hold on;
            p2=plot(sigma, perfMat(18:23,:));
            p3=plot(sigma, perfMat(24,:), 'Color', [0.9,0.5,0.2]);
            set(p1,'LineWidth',0.5);set(p2,'LineWidth',0.5);set(p3,'LineWidth',0.5);
            ylim([0.5,1]); box off
            fig4.cleanup; fig4.save([plots.dirs_save 'performance_cat4']);
%             legend('Opt','SumXT1','SumXT2', 'SumXT3','SumXT12','SumXT13','SumXT23', 'SumXT123');
%             xlabel('sigma');
%             ylabel('performance'); title('performance curve')
        end
        function prediction_relative_to_opt2_all
            model_idx = 1:25;
            sigma = 0;
            dirname_temp = compute.model_names{model_idx(1)};
            load([compute.dirs dirname_temp '/pred2opt/pred2opt.mat']);
            predMat = zeros(length(model_idx), length(pred2opt));
            for ii = 1:length(model_idx)
                dirname = compute.model_names{model_idx(ii)};
                load([compute.dirs dirname '/pred2opt/pred2opt.mat']);
                predMat(ii,:) = pred2opt;
            end
            figure_size = [40,28];
            fig1 = Figure(101,'size',figure_size);
            p1 = plot(sigma, predMat(1,:),'k');hold on;
            p2 = plot(sigma, predMat([2:5,25],:));
            set(p1,'LineWidth',.75); set(p2,'LineWidth',.75);
            ylim([0.5,1]); box off
            set(gca,'LineWidth',.5);
            
%             legend('Opt', 'Sum', 'Max', 'Min', 'Var');
%             xlabel('Sesory noise level');
%             ylabel('Proportion agreement with Opt'); 
%             title('prediction similarity to opt model, category 1');
            fig1.cleanup; fig1.save('pred2opt-cat1.eps');
            fig2 = Figure(102,'size',figure_size);
            plot(sigma, predMat(1,:), 'k'); hold on;
            plot(sigma, predMat(6:10,:));
            ylim([0.5,1]); box off
            set(gca,'LineWidth',.5);
%             legend('Opt','MaxT2','MaxT12','MaxT13','MaxT23','MaxT123');
%             xlabel('sigma');
%             ylabel('similarity to opt model'); title('prediction similarity to opt model, category 2');
            fig2.cleanup; fig2.save([plots.dirs_save 'pred2opt-cat2.eps']);
            fig3 = Figure(103,'size',figure_size);
            plot(sigma, predMat(1,:), 'k'); hold on;
            plot(sigma, predMat(11:16,:));
            plot(sigma, predMat(17,:), 'Color', [0.9,0.5,0.2]);
            ylim([0.5,1]); box off
            set(gca,'LineWidth',.5);
%             legend('Opt','SumErf','SumErfT1', 'SumErfT2','SumErfT3','SumErfT12','SumErfT13', 'SumErfT23');
%             xlabel('sigma');
%             ylabel('similarity to opt model'); title('prediction similarity to opt model, category 3');
            fig3.cleanup; fig3.save([plots.dirs_save 'pred2opt-cat3.eps']);
            fig4 = Figure(104,'size',figure_size);
            plot(sigma, predMat(1,:), 'k'); hold on;
            plot(sigma, predMat(18:23,:));
            plot(sigma, predMat(24,:), 'Color', [0.9,0.5,0.2]);
            ylim([0.5,1]); box off
            set(gca,'LineWidth',.5);
%             legend('Opt','SumXT1','SumXT2', 'SumXT3','SumXT12','SumXT13','SumXT23', 'SumXT123');
%             xlabel('sigma');
%             ylabel('similarity to opt model'); title('prediction similarity to opt model, category 3');
            fig4.cleanup; fig4.save([plots.dirs_save 'pred2opt-cat4.eps']);
        end
        function prediction_relative_to_opt2(model_idx)
            sigma = 0;
            dirname_temp = compute.model_names{model_idx(1)};
            load([compute.dirs dirname_temp '/pred2opt/pred2opt.mat']);
            predMat = zeros(length(model_idx), length(pred2opt));
            nameMat = cell(1,length(model_idx));
            for ii = 1:length(model_idx)
                dirname = compute.model_names{model_idx(ii)};
                load([compute.dirs dirname '/pred2opt/pred2opt.mat']);
                predMat(ii,:) = pred2opt;
                nameMat{ii} = dirname;
            end
            
            figure; plot(sigma, predMat(1,:),'k');hold on;
            plot(sigma, predMat(2:end,:));
            ylim([0.5,1]); box off
            legend(nameMat);
            xlabel('sigma');
            ylabel('similarity to opt model'); title('prediction similarity to opt model');
            
        end
        function performance_relative_to_opt2_all
            model_idx = 1:24;
            sigma = 0;
            dirname_temp = compute.model_names{model_idx(1)};
            load([compute.dirs dirname_temp '/pred2opt/pred2opt.mat']);
            predMat = zeros(length(model_idx), length(pred2opt));
            for ii = 1:length(model_idx)
                dirname = compute.model_names{model_idx(ii)};
                load([compute.dirs dirname '/pred2opt/pred2opt.mat']);
                predMat(ii,:) = pred2opt;
            end
            
            figure; plot(sigma, predMat(1,:),'k');hold on;
            plot(sigma, predMat(2:5,:));
            ylim([0.5,1]); box off
            legend('Opt', 'Sum', 'Max', 'Min', 'Var');
            xlabel('sigma');
            ylabel('similarity to opt model'); title('prediction similarity to opt model, category 1');
            figure;
            plot(sigma, predMat(1,:), 'k'); hold on;
            plot(sigma, predMat(6:10,:));
            ylim([0.5,1]); box off
            legend('Opt','Max2','MinMax','MinVar','MaxVar','MinMaxVar');
            xlabel('sigma');
            ylabel('similarity to opt model'); title('prediction similarity to opt model, category 2');
            figure;
            plot(sigma, predMat(1,:), 'k'); hold on;
            plot(sigma, predMat(11:16,:));
            plot(sigma, predMat(17,:), 'Color', [0.9,0.5,0.2]);
            ylim([0.5,1]); box off
            legend('Opt','SumErf','SumMin', 'SumMax','SumVar','SumMinMax','SumMinVar', 'SumMaxVar');
            xlabel('sigma');
            ylabel('similarity to opt model'); title('prediction similarity to opt model, category 3');
        end
        function model_similarity_mds
            nModels = 25;
            model_idx = 1:25;
            load('similarity/dissimilarity.mat');
            Y = mdscale(distMat,2);
            nameMat = cell(1,nModels);
            for ii = 1:nModels
                nameMat{ii} = plots.model_names{ii};
            end
            subj_idx = 1:9;
            bmcMat = zeros(1,nModels);
            for ii = 1:nModels
                dirname = compute.model_names{model_idx(ii)};
                eviMat = zeros(1,length(subj_idx));
                for jj = 1:length(subj_idx)
                    subjid = compute.subjids{subj_idx(jj)};
                    load([compute.dirs dirname '/evi/' subjid '.mat']);
                    eviMat(jj) = bmc;
                end
                bmcMat(ii) = mean(eviMat);
            end
            bmcMat = bmcMat - bmcMat(1);
            colormap_jet = colormap;
            % rescale likelihood vector on the colormap
%             bmcMat = (bmcMat - min(bmcMat))./(max(bmcMat) - min(bmcMat));
            idx = interp1(linspace(min(bmcMat), max(bmcMat),length(colormap_jet)),1:length(colormap_jet),bmcMat,'nearest','extrap');
            fig = Figure(101,'size',[150,80]);
            h=gscatter(Y(:,1),Y(:,2),1:25,colormap_jet(idx,:)); caxis([min(bmcMat), max(bmcMat)]);colorbar; hold on
            set(gca,'XTick',[],'YTick',[]); box on
            for ii = 1:nModels
                text(Y(ii,1),Y(ii,2),plots.model_names{ii});
            end
            fig.cleanup; fig.save('model_similarity_mds.eps');
            [idx,C] = kmeans(Y,3);
           
           
        end
        function model_similarity_matrix
            load('similarity/dissimilarity.mat');
            similarity = 1-distMat;
            figure; imagesc(similarity);caxis([0.5,1]);colorbar
            
            for ii = 1:24
                for jj = 1:24
                    text(ii-0.25,jj,num2str(similarity(ii,jj),2))
                end
            end
            nameMat = cell(1,24);
            for ii = 1:24
                nameMat{ii} = plots.model_names{ii};
            end
            set(gca, 'XTick',1:24)
            set(gca, 'XTickLabel',nameMat);
            
            set(gca, 'YTick',1:24)
            set(gca, 'YTickLabel',nameMat);
            set(gcf, 'Position', [1,1,1500,1000]);
            set(gca, 'TickLength',[0,0]);
            rotateXLabels(gca,90)
            set(gca,'XAxisLocation','bottom');
        end
        function compare_fit_pars(model_idx,tosave)
            % scatter plot to show the fit pars of each model vs optimal
            % model, each dot represents a subject. SS 2-27-2014
            if ~exist('tosave', 'var')
                tosave = 0;
            end
            for ii = model_idx
                dirname = compute.model_names{ii};
                optMat_lambda = zeros(1,length(compute.subjids));
                optMat_guess = zeros(1,length(compute.subjids));
                modelMat_lambda = zeros(1,length(compute.subjids));
                modelMat_guess = zeros(1,length(compute.subjids));
                for jj = 1:length(compute.subjids)
                    load([compute.dirs 'opt/fit_pars/' compute.subjids{jj}]);
                    optMat_lambda(jj) = CPG_lambda_hat;
                    optMat_guess(jj) = CPG_guess_hat;
                    load([compute.dirs dirname '/fit_pars/' compute.subjids{jj}]);
                    modelMat_lambda(jj) = CPG_lambda_hat;
                    modelMat_guess(jj) = CPG_guess_hat;
                                       
                end
                optMat_lambda
                modelMat_lambda
                figure; set(gcf, 'Position', get(gcf, 'Position').*[1,1,0.8,0.5]);
                subplot 121
                scatter(optMat_lambda, modelMat_lambda, 'o'); refline(1)
                title([dirname '-lambda']);
                subplot 122
                scatter(optMat_guess, modelMat_guess, 'o'); refline(1)
                title([dirname '-guess']);
                if tosave==1
                    save(a);
                end
            end
        end
        function bmc_agreement_opt(subj_idx)
            % plot mean bmc across subjects for a model as a function of
            % agreement to optimal model, compute R^2 of linear regression
            % SS 2014-03-04
            model_idx = [1:5,25,6:24];
            
            eviMat = zeros(length(subj_idx),length(model_idx));
            agreementMat = zeros(1,length(model_idx));
            for ii = 1:length(model_idx)
                dirname = compute.model_names{model_idx(ii)};
                
                for jj = 1:length(subj_idx)
                    subjid = compute.subjids{subj_idx(jj)};
                    load([compute.dirs dirname '/evi/' subjid '.mat']);
                    eviMat(jj,ii) = bmc;
                end
                
                load([compute.dirs dirname '/pred2opt/pred2opt.mat']);
                agreementMat(ii) = mean(pred2opt);
            end
            eviMat = bsxfun(@minus, eviMat, eviMat(:,1));
            bmcMat = mean(eviMat);
            steMat = std(eviMat)/sqrt(length(subj_idx));
            p = polyfit(agreementMat, bmcMat,1);
            yfit = polyval(p,agreementMat);
            yresid = bmcMat - yfit;
            SSresid = sum(yresid.^2);
            SStotal = (length(bmcMat)-1) * var(bmcMat);
            rsq = 1 - SSresid/SStotal;
            pearson_r = corrcoef(agreementMat, bmcMat);
            r = pearson_r(1,2)
            fig = Figure(101,'size',[70,50]); 
            h = errorbar(agreementMat(1),bmcMat(1),steMat(1),'LineStyle','None'); hold on
            set(h,'Color','k','Marker','o','MarkerSize',5);
            h = errorbar(agreementMat(2:6), bmcMat(2:6),steMat(2:6),'LineStyle','None');
            set(h,'Color','b','Marker','o','MarkerSize',5);
            h = errorbar(agreementMat(7:11), bmcMat(7:11),steMat(7:11),'LineStyle','None');
            set(h,'Color','r','Marker','o','MarkerSize',5);
            h = errorbar(agreementMat(12:25), bmcMat(12:25),steMat(12:25),'LineStyle','None');
            set(h,'Color','g','Marker','o','MarkerSize',5);
            
            legend('Opt', 'Class 1','Class 2','Class 3','Location','NorthWest');
            set(gca,'Clipping','off')
            [x,idx] = sort(agreementMat);
            plot(x, yfit(idx),'k--');
            xlim([0.6,1]); text(0.8,-1200,['r = ' num2str(r)]);
%             ylim([-1300,-900]);
%             xlabel('Model agreement relative to optimal model');
%             ylabel('Model log likelihood');
            fig.cleanup; fig.save([plots.dirs_save4 'bmc_agreement.eps']);
            
        end
        function bmc_agreement_all(subj_group_id,sz,tosave)
            
            % plot mean bmc across subjects for a model as a function of
            % agreement to each model, compute R^2 of linear regression
            % SS 2015-07-15
            if subj_group_id == 1 % task with correctness feedback
                subj_idx = 1:9;
            elseif subj_group_id == 2 % task without correctness feedback
                subj_idx = 10:14;
            else
                error('Invalid input for subject group id, please enter 1 or 2, 1 is with feedback, 2 is not:')
            end
            
            if ~exist('sz','var')
                sz = [200,175];
            end
            
            if ~exist('tosave','var')
                tosave = 0;
            end
                
            model_idx = [1:5,25,6:24];
            % load similarity mat
            load('similarity/dissimilarity.mat')
            similarity = similarity([1:5,25,6:24],[1:5,25,6:24]);
           
            % get agreement mat
            eviMat = utils.get_evi(subj_idx,model_idx,3);
            eviMat = eviMat/2000;
            bmcMat = mean(eviMat);
            steMat = std(eviMat)/sqrt(length(subj_idx));
            fig = Figure(140,'size',sz);
            for ii = 1:length(model_idx)
                subplot(5,5,ii)
                agreementMat = similarity(ii,:);
                p = polyfit(agreementMat, bmcMat,1);
                yfit = polyval(p,agreementMat);
                yresid = bmcMat - yfit;
                SSresid = sum(yresid.^2);
                SStotal = (length(bmcMat)-1) * var(bmcMat);
                rsq = 1 - SSresid/SStotal;
                pearson_r = corrcoef(agreementMat, bmcMat);
                r = pearson_r(1,2);
                % save rsq as global maximum index
                dirname = compute.model_names{model_idx(ii)};
                save_dir = ['similarity/' dirname '_rsq_real_opt_group' num2str(subj_group_id)];
                
                save(save_dir, 'rsq','r');
                
                h = errorbar(agreementMat(1),bmcMat(1),steMat(1),'LineStyle','None'); hold on
                set(h,'Color','k','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(2:6), bmcMat(2:6),steMat(2:6),'LineStyle','None');
                set(h,'Color','b','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(7:11), bmcMat(7:11),steMat(7:11),'LineStyle','None');
                set(h,'Color','r','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(12:25), bmcMat(12:25),steMat(12:25),'LineStyle','None');
                set(h,'Color','g','Marker','o','MarkerSize',2.5);
%                 if ii==1
%                     legend('Opt', 'Class 1','Class 2','Class 3','Location','NorthWest');
%                 end
                set(gca,'Clipping','off')
                [x,idx] = sort(agreementMat);
                plot(x, yfit(idx),'k--');
                ylim([-0.7,-0.3])
                xlim([0.6,1]); text(0.7,-0.4,['\it r = ' num2str(r,'%.02f')]);
                title(plots.model_names{model_idx(ii)})
            end
            fig.cleanup;
            if tosave
                fig.save([plots.dirs_save4 'bmc_agreement_all_real_group' num2str(subj_group_id) '.eps']);
            end
            
        end
        function bmc_agreement_fake_all(runs,sz,tosave)
            % plot mean bmc across fake subjects as a function of
            % agreement to optimal model, compute R^2 of linear regression
            % SS 2015-07-08
            model_idx = [1:5,25,6:24];
            if ~exist('sz','var');
                sz = [200,175];
            end
            if ~exist('tosave','var')
                tosave = 0;
            end
            fig = Figure(120,'size',sz);
            for kk = 1:length(model_idx)
                subplot(5,5,kk)
                subj_ids = utils.create_subjids(model_idx(kk),runs);
                eviMat = zeros(length(subj_ids),length(model_idx));
                agreementMat = zeros(1,length(model_idx));
                for ii = 1:length(model_idx)
                    dirname = compute.model_names{model_idx(ii)};

                    for jj = 1:length(subj_ids)
                        subjid = subj_ids{jj};
                        load([compute.dirs_fake1 dirname '/evi/' subjid '.mat']);
                        eviMat(jj,ii) = bmc;
                    end
                    
                    load([compute.dirs dirname '/pred2opt/pred2opt.mat']);
                    agreementMat(ii) = mean(pred2opt);
                end
                eviMat = eviMat/3000;
                bmcMat = mean(eviMat);
                steMat = std(eviMat)/sqrt(length(subj_ids));
                p = polyfit(agreementMat, bmcMat,1);
                yfit = polyval(p,agreementMat);
                yresid = bmcMat - yfit;
                SSresid = sum(yresid.^2);
                SStotal = (length(bmcMat)-1) * var(bmcMat);
                rsq = 1 - SSresid/SStotal;
                pearson_r = corrcoef(agreementMat, bmcMat);
                r = pearson_r(1,2);
                % save rsq as global maximum index
                dirname = compute.model_names{model_idx(kk)};
                save_dir = ['similarity/' dirname '_rsq_fake_opt.mat'];
                
                save(save_dir, 'rsq','r');
                
                h = errorbar(agreementMat(1),bmcMat(1),steMat(1),'LineStyle','None'); hold on
                set(h,'Color','k','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(2:6), bmcMat(2:6),steMat(2:6),'LineStyle','None');
                set(h,'Color','b','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(7:11), bmcMat(7:11),steMat(7:11),'LineStyle','None');
                set(h,'Color','r','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(12:25), bmcMat(12:25),steMat(12:25),'LineStyle','None');
                set(h,'Color','g','Marker','o','MarkerSize',2.5);
%                 if kk==1
%                     legend('Opt', 'Class 1','Class 2','Class 3','Location','NorthWest');
%                 end
                set(gca,'Clipping','off')
                [x,idx] = sort(agreementMat);
                plot(x, yfit(idx),'k--');
%                 yLim = [min(yfit)-50,max(yfit)+50];
                ylim([-0.7,-0.3])
                xlim([0.6,1]); text(0.75,-0.35,['\it r\rm = ' num2str(r,'%.02f')]);
                title(plots.model_names{model_idx(kk)})
            end
            fig.cleanup
            if tosave
                fig.save([plots.dirs_save4 'bmc_agreement_fake_opt.eps']);
            end
        end
        function bmc_agreement_match(runs,sz,tosave)
            model_idx = [1:5,25,6:24];
            if ~exist('sz','var');
                sz = [200,175];
            end
            if ~exist('tosave','var')
                tosave = 0;
            end
            load('similarity/dissimilarity.mat')
            similarity = similarity([1:5,25,6:24],[1:5,25,6:24]);
            eviMat_all = utils.get_evi(model_idx,model_idx,3,0,runs);
            fig = Figure(150,'size',sz);
            for kk = 1:length(model_idx)
                subplot(5,5,kk)
                eviMat = squeeze(eviMat_all(kk,:,:));
                eviMat = eviMat/3000;
                agreementMat = similarity(kk,:);
                bmcMat = mean(eviMat);
                steMat = std(eviMat)/sqrt(length(runs));
                p = polyfit(agreementMat, bmcMat,1);
                yfit = polyval(p,agreementMat);
                yresid = bmcMat - yfit;
                SSresid = sum(yresid.^2);
                SStotal = (length(bmcMat)-1) * var(bmcMat);
                rsq = 1 - SSresid/SStotal;
                pearson_r = corrcoef(agreementMat, bmcMat);
                r = pearson_r(1,2);
                % save rsq as global maximum index
                dirname = compute.model_names{model_idx(kk)};
                save_dir = ['similarity/' dirname '_rsq_fake_match'];
                
                save(save_dir, 'rsq','r');
                
                h = errorbar(agreementMat(1),bmcMat(1),steMat(1),'LineStyle','None'); hold on
                set(h,'Color','k','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(2:6), bmcMat(2:6),steMat(2:6),'LineStyle','None');
                set(h,'Color','b','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(7:11), bmcMat(7:11),steMat(7:11),'LineStyle','None');
                set(h,'Color','r','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(12:25), bmcMat(12:25),steMat(12:25),'LineStyle','None');
                set(h,'Color','g','Marker','o','MarkerSize',2.5);
%                 if kk==1
%                     legend('Opt', 'Class 1','Class 2','Class 3','Location','NorthWest');
%                 end
                set(gca,'Clipping','off')
                [x,idx] = sort(agreementMat);
                plot(x, yfit(idx),'k--');
                ylim([-0.7,-0.3])
                xlim([0.6,1]); text(0.65,-0.35,['\it r\rm = ' num2str(r,'%.02f')]);
                title(plots.model_names{model_idx(kk)})
            end
            fig.cleanup
            if tosave
                fig.save([plots.dirs_save4 'bmc_agreement_fake_match.eps']);
            end
        end
        function bmc_agreement_summary_real(subj_group_id,tosave)
            
            if ~exist('tosave','var')
                tosave = 0;
            end
            similar_idx = [1,14,16,17,20,22,23,24];
            idx_all = 1:25;
            non_similar_idx = idx_all(~ismember(idx_all,similar_idx));
            
            rsqMat = zeros(size(idx_all));
            for ii = 1:length(idx_all)
                dirname = compute.model_names{idx_all(ii)};
                load(['similarity/' dirname '_rsq_real_opt_group' num2str(subj_group_id) '.mat'])
                rsqMat(ii) = r;
            end
            
            similar_rsq = rsqMat(similar_idx);
            non_similar_rsq = rsqMat(non_similar_idx);
            
            p = ranksum(similar_rsq, non_similar_rsq)
            fig = Figure(200,'size',[50,40]); hold on
            
            plot([0.85,1.15],[similar_rsq(1),similar_rsq(1)],'r');
            plot(1,similar_rsq(2:end),'ko');
            for ii = 1:length(non_similar_rsq)
                plot([1.9,2.1],[non_similar_rsq(ii), non_similar_rsq(ii)], 'k');
            end
            set(gca, 'XTick', [1,2]);
            set(gca, 'XTickLabel',{'similar', 'dissimilar'})
            xlim([0,3])
            fig.cleanup
            
            if tosave
                fig.save([plots.dirs_save4 'gmi_real_group_' num2str(subj_group_id) '.eps'])
            end

            
        end
        function bmc_agreement_summary_fake_opt(tosave)
            
            if ~exist('tosave','var')
                tosave = 0;
            end
            similar_idx = [1,14,16,17,20,22,23,24];
            idx_all = 1:25;
            non_similar_idx = idx_all(~ismember(idx_all,similar_idx));
            
            rsqMat = zeros(size(idx_all));
            for ii = 1:length(idx_all)
                dirname = compute.model_names{idx_all(ii)};
                load(['similarity/' dirname '_rsq_fake_opt.mat'])
                rsqMat(ii) = r;
            end
           
            similar_rsq = rsqMat(similar_idx);
            non_similar_rsq = rsqMat(non_similar_idx);
            p = ranksum(similar_rsq, non_similar_rsq);
            fig = Figure(200,'size',[50,40]); hold on
            plot(1,similar_rsq,'ko');
            for ii = 1:length(non_similar_rsq)
                plot([1.9,2.1],[non_similar_rsq(ii),non_similar_rsq(ii)],'k');
            end
            set(gca, 'XTick', [1,2]);
            set(gca, 'XTickLabel',{'similar', 'dissimilar'})
            xlim([0,3])
            ylim([-1,1])
            fig.cleanup
            if tosave
                fig.save([plots.dirs_save4 'gmi_fake.eps'])
            end
        end
        function bmc_agreement_summary_fake_match(tosave)
            if ~exist('tosave','var')
                tosave = 0;
            end
            idx_all = 1:25;
            non_similar_idx = [2:5,25,6:13,15,18:19,21];
            
            rsqMat1 = zeros(size(idx_all));
            rsqMat2 = zeros(size(idx_all));
            for ii = 1:length(idx_all)
                dirname = compute.model_names{idx_all(ii)};
                load(['similarity/' dirname '_rsq_real_opt.mat']);
                rsqMat1(ii) = r;
                load(['similarity/' dirname '_rsq_fake_match.mat'])
                rsqMat2(ii) = r;
            end
           
            models = plots.model_names(non_similar_idx);
           
            non_similar_rsq_real = rsqMat1(non_similar_idx);
            non_similar_rsq_fake = rsqMat2(non_similar_idx);
            p = signrank(non_similar_rsq_real, non_similar_rsq_fake);
            fig = Figure(200,'size',[150,40]); hold on
            for ii = 1:length(non_similar_rsq_real)
                plot([ii-0.15,ii+.15], [non_similar_rsq_real(ii),non_similar_rsq_real(ii)],'k');
            end
            plot(1:length(non_similar_rsq_fake), non_similar_rsq_fake,'ko');
            set(gca, 'xTick',1:length(non_similar_rsq_real))
            set(gca, 'xTickLabel',models)
            
            xlim([0,length(non_similar_idx)+2])
            fig.cleanup
            if tosave
                fig.save([plots.dirs_save4 'gmi_match.eps'])
            end
        end
        function gof_mll(subj_idx)
            % relationship between goodness of fit and model log likelihood
            % SS 2015-07-22
            model_idx = 1:25;
            subjid = compute.subjids{subj_idx};
            [~,response] = utils.readdata(subjid);
            predictionMat = zeros(2000,length(model_idx));
            eviMat = zeros(1,length(model_idx));
            % load fit results
            for ii = model_idx
                dirname = compute.model_names{ii};
                load([compute.dirs dirname '/fit_results/' subjid '.mat']);
                predictionMat(:,ii) = CPG_prediction;
                load([compute.dirs dirname '/evi/' subjid '.mat']);
                eviMat(ii) = bmc;
            end
            predictionMat(response==-1,:) = 1-predictionMat(response==-1,:);
            gof = mean(predictionMat);
            figure; plot(gof, eviMat,'o')
        end
        function gof_agreement_all(subj_idx,sz,tosave)
            % plot mean gof across subjects for a model as a function of
            % agreement to each model, compute R^2 of linear regression
            % SS 2015-07-15
            
            if ~exist('sz','var')
                sz = [200,175];
            end
            
            if ~exist('tosave','var')
                tosave = 0;
            end
                
            model_idx = [1:5,25,6:24];
            % load similarity mat
            load('similarity/dissimilarity.mat')
            similarity = similarity([1:5,25,6:24],[1:5,25,6:24]);
           
            % get prediction mat
            predictionMat = zeros(length(subj_idx),length(model_idx));
            
            for ii = 1:length(subj_idx)
                subjid = compute.subjids{subj_idx(ii)};
                [~,response] = utils.readdata(subjid);
                for jj = 1:length(model_idx)
                    dirname = compute.model_names{model_idx(jj)};
                    load([compute.dirs dirname '/fit_results/' subjid '.mat']);
                    CPG_prediction(response==-1) = 1 - CPG_prediction(response==-1);
                    predictionMat(ii,jj) = mean(CPG_prediction);
                    
                end
            end
            
            gof = mean(predictionMat);
            ste_gof = std(predictionMat)/sqrt(length(subj_idx));
            fig = Figure(140,'size',sz);
            for ii = 1:length(model_idx)
                subplot(5,5,ii)
                agreementMat = similarity(ii,:);
                p = polyfit(agreementMat, gof,1);
                yfit = polyval(p,agreementMat);
                yresid = gof - yfit;
                SSresid = sum(yresid.^2);
                SStotal = (length(gof)-1) * var(gof);
                rsq = 1 - SSresid/SStotal;
                
                h = errorbar(agreementMat(1),gof(1),ste_gof(1),'LineStyle','None'); hold on
                set(h,'Color','k','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(2:6), gof(2:6),ste_gof(2:6),'LineStyle','None');
                set(h,'Color','b','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(7:11), gof(7:11),ste_gof(7:11),'LineStyle','None');
                set(h,'Color','r','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(12:25), gof(12:25),ste_gof(12:25),'LineStyle','None');
                set(h,'Color','g','Marker','o','MarkerSize',2.5);
%                 if ii==1
%                     legend('Opt', 'Class 1','Class 2','Class 3','Location','NorthWest');
%                 end
                set(gca,'Clipping','off')
                [x,idx] = sort(agreementMat);
                plot(x, yfit(idx),'k--');
                ylim([0.5,0.8])
                xlim([0.6,1]); text(0.7,0.75,['\it R\rm^{2} = ' num2str(rsq,'%.02f')]);
                title(plots.model_names{model_idx(ii)})
            end
            fig.cleanup;
            if tosave
                fig.save([plots.dirs_save4 'gof_agreement_all_real.eps']);
            end
        end
        function gof_agreement_match(runs,sz,tosave)
            model_idx = [1:5,25,6:24];
            if ~exist('sz','var');
                sz = [200,175];
            end
            if ~exist('tosave','var')
                tosave = 0;
            end
            load('similarity/dissimilarity.mat')
            similarity = similarity([1:5,25,6:24],[1:5,25,6:24]);
            % get prediction mat
            predictionMat = zeros(length(model_idx),length(runs),length(model_idx));
            
            for ii = 1:length(model_idx) % different data
                name_base = compute.model_names{model_idx(ii)};
                for jj = 1:length(runs)
                     subjid = [name_base '_' num2str(runs(jj))];
                     [~,response] = utils.readfakedata(subjid);
                    for kk = 1:length(model_idx) % different models
                        dirname = compute.model_names{model_idx(kk)};
                        filename = [compute.dirs_fake1 dirname '/fit_results/' subjid '.mat'];
                        if exist(filename, 'file')
                            load([compute.dirs_fake1 dirname '/fit_results/' subjid '.mat']);
                            CPG_prediction(response==-1) = 1-CPG_prediction(response==-1);
                            predictionMat(ii,jj,kk) = mean(CPG_prediction);
                        else
                        predictionMat(ii,jj,kk) = 0;
                        end
                    end
                end
            end
            
           
            fig = Figure(150,'size',sz);
            for kk = 1:length(model_idx)
                subplot(5,5,kk)
                
                agreementMat = similarity(kk,:);
                gofMat = squeeze(predictionMat(kk,:,:));
                gof = mean(gofMat);
                ste_gof = std(gofMat)/sqrt(size(gofMat,1));
                p = polyfit(agreementMat, gof,1);
                yfit = polyval(p,agreementMat);
                yresid = gof - yfit;
                SSresid = sum(yresid.^2);
                SStotal = (length(gof)-1) * var(gof);
                rsq = 1 - SSresid/SStotal;
                   
                h = errorbar(agreementMat(1),gof(1),ste_gof(1),'LineStyle','None'); hold on
                set(h,'Color','k','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(2:6), gof(2:6),ste_gof(2:6),'LineStyle','None');
                set(h,'Color','b','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(7:11), gof(7:11),ste_gof(7:11),'LineStyle','None');
                set(h,'Color','r','Marker','o','MarkerSize',2.5);
                h = errorbar(agreementMat(12:25), gof(12:25),ste_gof(12:25),'LineStyle','None');
                set(h,'Color','g','Marker','o','MarkerSize',2.5);
%                 if kk==1
%                     legend('Opt', 'Class 1','Class 2','Class 3','Location','NorthWest');
%                 end
                set(gca,'Clipping','off')
                [x,idx] = sort(agreementMat);
                plot(x, yfit(idx),'k--');
                ylim([0.5,0.8])
                xlim([0.6,1]); text(0.65,0.75,['\it R\rm^{2} = ' num2str(rsq,'%.02f')]);
                title(plots.model_names{model_idx(kk)})
            end
            fig.cleanup
            if tosave
                fig.save([plots.dirs_save4 'gof_agreement_fake_match.eps']);
            end
        end
        function simulation_learning
            sigma = 0;
            load('learning.mat'); nTrials = length(sigEstMat);
            sigEstMean = mean(sigEstMat); sigEstStd = std(sigEstMat);
            sigsEstMean = mean(sigsEstMat); sigsEstStd = std(sigsEstMat);
            figure; set(gcf, 'Position', get(gcf,'Position').*[1,1,1,0.5]);
            subplot 121
            plot(1:nTrials,sigEstMean); hold on
            plot(1:nTrials,sigEstMean+sigEstStd,'b:');
            plot(1:nTrials,sigEstMean-sigEstStd,'b:');
            plot(1:nTrials,sigma*ones(1,nTrials),'r');
            legend('mean','mean+sem','mean-sem','true value');
            title('Estimated sigma across trials');
            xlabel('Trial number'); ylabel('Estimated sigma');
            subplot 122
            plot(log(1:nTrials),sigsEstMean); hold on
            plot(log(1:nTrials),sigsEstMean+sigsEstStd,'b:');
            plot(log(1:nTrials),sigsEstMean-sigsEstStd,'b:');
            plot(log(1:nTrials),compute.s_std*ones(1,nTrials),'r');
            legend('mean','mean+std','mean-std','true value');
            title('Estimated sigma s across trials');
            xlabel('Trial number'); ylabel('Estimated sigma s');
        end
        function performance_trials(subj_idx,window_size)
            % load data
            subjid = compute.subjids{subj_idx};
            files = dir(['session1/output/' subjid '*.mat']);
            % read every file
            stimulus = [];
            response = [];
            performance = [];
            for ii=1:length(files)
                load(['session1/output/' files(ii).name]);
                stimulus = [stimulus; data(:,1:2)];
                response = [response; data(:,3)];
                performance = [performance; data(:,4)];
            end
            perf = zeros(size(window_size+1:length(performance)-window_size));
            cnt = 0;
            for ii = window_size+1:length(performance)-window_size
                cnt = cnt+1;
                perf(cnt) = mean(performance(ii-window_size:ii+window_size));
            end
            trial_no = window_size+1:length(performance)-window_size;
            figure; plot(trial_no,perf); ylim([0.5,1]); xlim([min(trial_no), max(trial_no)]);
        end
        function performance_trials_all(window_size,session)
            % load data
            file = ['output2/session' num2str(session) '/' compute.subjids{1} '_' num2str(session) '.mat'];
            load(file); performance = data(:,4);
            perfMat = zeros(length(compute.subjids),length(window_size+1:length(performance)-window_size));
            for jj = 1:length(compute.subjids)
                subjid = compute.subjids{jj};
                files = dir(['output2/session' num2str(session) '/' subjid '*.mat']);
                % read every file
                stimulus = [];
                response = [];
                performance = [];
                for ii=1:length(files)
                    load(['output2/session' num2str(session) '/' files(ii).name]);
                    stimulus = [stimulus; data(:,1:2)];
                    response = [response; data(:,3)];
                    performance = [performance; data(:,4)];
                end
                perf = zeros(size(window_size+1:length(performance)-window_size));
                cnt = 0;
                for ii = window_size+1:length(performance)-window_size
                    cnt = cnt+1;
                    perf(cnt) = mean(performance(ii-window_size:ii+window_size));
                end
                perfMat(jj,:)=perf;
            end
            trial_no = window_size+1:length(performance)-window_size;
            perf_mean = mean(perfMat);
            perf_sem = std(perfMat)/sqrt(length(compute.subjids));
            perf_up = perf_mean + perf_sem;
            perf_low = perf_mean - perf_sem;
            figure;  color = get(gca,'ColorOrder'); color = min(color(1,:)+ 0.65,1);
            patch([trial_no wrev(trial_no)], [perf_up wrev(perf_low)], color ,'LineStyle','None');
            hold on; plot(trial_no,perf_mean);
            ylim([0.5,1]); xlim([min(trial_no), max(trial_no)]); 
            title(['performance-session' num2str(session)]);
            xlabel('Number of trial'); ylabel('Performance');
        end
        function learning_curve(subj_idx,window_size)
            % load data
            subjid = compute.subjids{subj_idx};
            common_dir = compute.dirs3; dirname = 'opt';
            load([common_dir dirname '/fit_results_learning/' subjid]);
            perf = zeros(size(window_size+1:length(performance)-window_size));
            perf_pred = zeros(size(perf));
            cnt = 0;
            for ii = window_size+1:length(performance)-window_size
                cnt = cnt+1;
                perf(cnt) = mean(performance(ii-window_size:ii+window_size));
                perf_pred(cnt) = mean(CPG_prediction(ii-window_size:ii+window_size));
            end
            trial_no = window_size+1:length(performance)-window_size;
            figure; plot(trial_no,perf, trial_no, perf_pred); ylim([0.5,1]); xlim([min(trial_no), max(trial_no)]);
        end
        function learning_curve_all(window_size,session)
            dirname = compute.model_names{1};
            nSubjs = length(compute.subjids);
            common_dir = compute.dirs3;
            load([common_dir dirname '/fit_results_learning/' compute.subjids{1} '_' num2str(session)]);
            trial_no = window_size+1:length(performance)-window_size;
            perfMat = zeros(nSubjs,length(trial_no));
            perf_predMat = zeros(nSubjs,length(trial_no));
            for jj = 1:nSubjs
                cnt = 0;
                load([common_dir dirname '/fit_results_learning/' compute.subjids{jj} '_' num2str(session)]);
                for ii = window_size+1:length(performance)-window_size
                    cnt = cnt+1;
                    perfMat(jj,cnt) = mean(performance(ii-window_size:ii+window_size));
                    perf_predMat(jj,cnt) = mean(CPG_prediction(ii-window_size:ii+window_size));
                end
            end
            perf = squeeze(mean(perfMat));
            perf_pred = squeeze(mean(perf_predMat));
            err_data = squeeze(std(perfMat)/sqrt(nSubjs));
            err_pred     = squeeze(std(perf_predMat,1)/sqrt(nSubjs));
            upper_pred   = perf_pred + err_pred;
            lower_pred   = perf_pred - err_pred;
            upper_data = perf+err_data;
            lower_data = perf-err_data;


            a = figure; 
            set(gcf, 'Position', get(gcf, 'Position').*[1,1,0.8,0.8]);

            colorvec = get(gca, 'ColorOrder');
            colorvec = min(colorvec+.65,1);
            patch([trial_no wrev(trial_no)],[upper_pred wrev(lower_pred)],  colorvec(7,:), 'Linestyle','None'); hold on;
            patch([trial_no wrev(trial_no)],[upper_data wrev(lower_data)],  colorvec(1,:), 'Linestyle','None');
            legend('Prediction','Data'); ylim([0.5,1])
            xlabel('Trial number'); ylabel('Performance');
            title(['learning curve, session' num2str(session)]);
        end
        function pars_regression(model_idx, run_idx, tosave)
            % parameter regression for fake data test
            if ~exist('tosave','var')
                tosave = 0;
            end
            for ii = model_idx
                dirname = compute.model_names{ii};
                % get real parameters and fit parameters
                real_pars = zeros(2,length(run_idx));
                fit_pars = zeros(2,length(run_idx));
                for jj = run_idx
                    subjid = [compute.model_names{ii} '_' num2str(jj)];
                    load([compute.dirs_fake1 dirname '/fit_pars/' subjid]);
                    fit_pars(1,jj) = CPG_lambda_hat;
                    fit_pars(2,jj) = CPG_guess_hat;
                    load([compute.dirs_fakedata1 subjid], 'lambda','guess_rate');
                    real_pars(1,jj) = lambda;
                    real_pars(2,jj) = guess_rate;
                end
                a = figure; set(gcf,'Position',[100,100,400,160]); 
                subplot 121; plot(real_pars(1,:), fit_pars(1,:),'k.','MarkerSize',15); xlim([0,0.3]);ylim([0,0.3]); hold on; h = refline(1); box off
                set(h,'Color','r','LineWidth',1);
                x = real_pars(1,:); y = fit_pars(1,:);
                y_pred = x;
                R21 = 1-sum((y-y_pred).^2)/sum((y-mean(y)).^2)
                set(gca, 'XTickLabel',[]);
                set(gca, 'YTickLabel',[]);
%                 xlabel('real lambda'); ylabel('fit lambda');
                subplot 122; plot(real_pars(2,:), fit_pars(2,:),'k.','MarkerSize',15); hold on; h=refline(1); box off
                set(h,'Color','r','LineWidth',1);
                x = real_pars(2,:); y = fit_pars(2,:);
                y_pred = x;
                R22 = 1-sum((y-y_pred).^2)/sum((y-mean(y)).^2)
                set(gca, 'XTickLabel',[]);
                set(gca, 'YTickLabel',[]);
%                 xlabel('real guessing rate'); ylabel('fit guessing rate');
                title(compute.model_names{ii})
                if tosave
                    saveas(a, [plots.dir_save 'par_regression_' dirname], 'emf');
                end
            end
        end
        function show_decision_rule_3d(model_idx)
            x1 = -20:0.5:20;
            x2 = -20:0.5:20;
            x3 = -20:0.5:20;
            nTrials = length(x1)*length(x2)*length(x3);
            c1 = ones(1,nTrials);
            c2 = ones(nTrials,3);
            xMat = ones(3,length(x1),length(x2),length(x3));
            a = 0.5*ones(length(x1),length(x2),length(x3));
            for ii = 1:length(x1)
                for jj = 1:length(x2)
                    for kk = 1:length(x3)
                        xMat(:,ii,jj,kk)=[x1(ii),x2(jj),x3(kk)];                        
                    end
                end
            end
            
            xMat2 = reshape(xMat,3,nTrials);
            obs_response = compute.decision_rule(model_idx,xMat2,0.04,nTrials);
            
            c1(obs_response>0) = 0;
            c2(obs_response>0,:) = repmat([0.6,0.6,1],sum(obs_response>0),1);
            
            c1 = reshape(c1,length(x1),length(x2),length(x3));
            c2 = reshape(c2,length(x1),length(x2),length(x3),3);
            save('opt_structure','c1','c2','x1','x2','x3');
            figure; 
            vol3d('xdata',x1,'ydata',x2,'zdata',x3,'cdata',c2); hold on
            p = patch(isosurface(x1,x2,x3,c1,0),'FaceColor','interp','EdgeColor','interp');
            isonormals(x1,x2,x3,c1,p)
            p.FaceColor = 'blue';
            p.EdgeColor = 'none';
            camlight
            lighting gouraud
            view(3)
            xlabel('x1'); ylabel('x2'); zlabel('x3')
        end
        function entropy_LML(subj_idx,model_idx)
            entropyMat_subj = zeros(1,length(subj_idx));
            LMLMat = zeros(length(subj_idx),length(model_idx));
            
            models = plots.model_names(model_idx);
            
            for ii = 1:length(subj_idx)
                subjid = compute.subjids{subj_idx(ii)};
                load(['real_data_results/entropy/' subjid '_9.mat']);
                entropyMat_subj(ii) = entropyvalue;
                for jj = 1:length(model_idx)
                    dirname = compute.model_names{model_idx(jj)};
                    load([compute.dirs dirname '/evi_bin/' subjid])
                    LMLMat(ii,jj) = bmc;
                end
            end
            LML_mean = mean(LMLMat);
            LML_ste = std(LMLMat)/sqrt(size(LMLMat,1));
            entropy_mean = mean(entropyMat_subj);
            entropy_ste = std(entropyMat_subj)/sqrt(size(LMLMat,1));
            [LML_sort_mean, sort_idx] = sort(LML_mean);
            LML_sort_mean = wrev(LML_sort_mean);
            sort_idx = wrev(sort_idx);
            LML_ste = LML_ste(sort_idx);
            models = models(sort_idx);
            models = {'Entropy',models{:}};
            fig = Figure(180,'size',[150,50]);
            bar(1:length(model_idx)+1,[entropy_mean,LML_sort_mean],'FaceColor',[0.7,0.7,0.7]); hold on
            errorbar(1:length(model_idx)+1,[entropy_mean,LML_sort_mean],[entropy_ste, LML_ste],'k','MarkerSize',2,'LineStyle','None');
           
            xlim([0,length(model_idx)+2])
            ylim([0.5,1]);
            set(gca, 'XTick', 1:length(model_idx)+1)
            set(gca, 'xTickLabel', models)
            
            fig.cleanup
            fig.save('Prediction performance')
            
        end
        function entropy_bins(subj_idx)
            nBins = [3,5,7,9,11];
            entropyMat_subj = zeros(length(subj_idx),length(nBins));
            for ii = 1:length(subj_idx)
                subjid = compute.subjids{ii};
                for jj = 1:length(nBins)
                    nbin = nBins(jj);
                    load(['real_data_results/entropy/' subjid '_' num2str(nbin)])
                    entropyMat_subj(ii,jj) = entropyvalue;
                end
            end
            entropyMat_subj = exp(entropyMat_subj/2000);
            figure; plot(nBins, entropyMat_subj, 'o-')
            xlabel('number of bins')
            ylabel('predictability')
            set(gca, 'XTick', nBins)
            xlim([2,12]); ylim([0.5,0.8])
        end
        function entropy_mle_subj(subj_idx,model_idx,cross)
            % For individual subjects, show mle and entropy of each model as a
            % horizontal bar
            % cross specifies whether cross validation is applied
            % modified 2015-09-14 SS
            if ~exist('cross','var')
                cross = 0;
            end
            
            if cross == 0
                mle_dir = '/mle_bin_grid/';
                entropy_dir = 'entropy_est_0/';
                nTrials = 2000;
                lower = -1500;
            else
                mle_dir = '/mle_bin_grid_cross/';
                entropy_dir = 'entropy_est_2/';
                nTrials = 1000;
                lower = -800;
            end
                
            fig = Figure(700,'size',[10*length(subj_idx)+10,60]); hold on
            chance = log(0.5)*nTrials;
            pred_prob = zeros(1,length(subj_idx));
            opt_prob = zeros(1,length(subj_idx));
            lml = zeros(1,length(subj_idx));
            err = zeros(1,length(subj_idx));
            Dkl_table = zeros(length(subj_idx),length(model_idx));
            entropy_table = zeros(length(subj_idx),1);
            mle_table = zeros(length(subj_idx),length(model_idx));
            p_table = zeros(length(subj_idx),length(model_idx));
            for ii = 1:length(subj_idx)
                subjid = compute.subjids{subj_idx(ii)};
                for jj = 1:length(model_idx)
                    dirname = compute.model_names{model_idx(jj)};
                    load([compute.dirs dirname mle_dir subjid])
                    Dkl_table(ii,jj) = dkl;
                    mle_table(ii,jj) = mle_mean;
                    p_table(ii,jj) = p;
                    if model_idx(jj)~=1 % suboptimal models
                        plot([ii-0.2,ii+0.2],[mle_mean,mle_mean],'color',[0.55,0.55,0.55],'linewidth',0.5)
                    else % optimal model
                        lml(ii) = mle_mean;
                        err(ii) = mle_err;
                        opt_prob(ii) = exp(mle_mean/nTrials);
                        exp(mle_range/nTrials)
                    end
                end
                plot([ii-0.3,ii+0.3],[lml(ii),lml(ii)],'b:','linewidth',0.5)
                
                load(['real_data_results/' entropy_dir subjid '_9.mat'])
                plot([ii-0.3,ii+0.3],[entropyvalue,entropyvalue],'g--','linewidth',0.5)
                
                pred_prob(ii) = exp(entropyvalue/nTrials);
                entropy_table(ii) = entropyvalue;
                
            end
            pred_prob
            errorbar(1:length(subj_idx), lml, err, 'b','LineStyle','None')
            ylim([lower,0])
            set(gca, 'xtick', 1:length(subj_idx),'XAxisLocation','top')
            plot(0.8:0.1:length(subj_idx)+0.2,chance*ones(1,length(0.8:0.1:length(subj_idx)+0.2)),'k.','linewidth',0.5)
            
            xlabel('Subject id'); ylabel('Negative entropy or LML')
            mean(pred_prob)
            std(pred_prob)/sqrt(length(subj_idx))
            
            mean(opt_prob)
            std(opt_prob)/sqrt(length(subj_idx))
            fig.cleanup
            fig.save('Entropy and LLmax.eps')
            pMat = p_table;
            p_table(length(subj_idx)+1,:) = mean(pMat);
            p_table(length(subj_idx)+2,:) = std(pMat)/sqrt(length(subj_idx));
%             T = array2table(round(p_table,2)', 'RowNames', plots.model_names(model_idx), 'VariableNames', [compute.subjids(subj_idx),'mean','sem'])
% 
%             writetable(T,'p_values.csv')
            pVec = zeros(1,length(model_idx));
            for ii = 1:length(model_idx)
                pVec(ii) = signrank(mle_table(:,ii),entropy_table,'tail','left');
            end
            pVec
            fig2 = Figure(800, 'size', [150,40]); hold on
            [sorted_pVec, idx] = sort(pVec,'descend');
            models = plots.model_names(model_idx(idx));
            plot(1:length(models), sorted_pVec, 'ko');
            set(gca, 'xTick',1:length(model_idx))
            set(gca, 'xTickLabel',models)
            fig2.cleanup
        end
        function entropy_sigma(subj_idx)
            % scatter plot of estimated entropy v.s. sigma across subjects
            est_neg_entropy = zeros(1,length(subj_idx));
            est_sigma = zeros(1,length(subj_idx));
            est_guess = zeros(1,length(subj_idx));
            for ii = 1:length(subj_idx)
                subjid = compute.subjids{subj_idx(ii)};
                dirname = 'opt';
                load([compute.dirs dirname '/mle_bin_grid_cross/' subjid '.mat'],'lambda_hat','guess_hat');
                load(['real_data_results/entropy_est_2/' subjid '_9.mat'],'entropyvalue');
                est_sigma(ii) = 1/sqrt(lambda_hat);
                est_neg_entropy(ii) = entropyvalue;
                est_guess(ii) = guess_hat;
                
            end
            fig = Figure(130,'size',[100,40]);
            subplot 121
            scatter(est_sigma, est_neg_entropy, 'ko')
            xlabel('estimated\sigma')
            ylabel('estimated negative entropy')
            subplot 122
            scatter(est_guess, est_neg_entropy, 'ko')
            xlabel('estimated lapse rate')
            fig.cleanup
            fig.save('entropy_sigma.eps')
        end
    end
end

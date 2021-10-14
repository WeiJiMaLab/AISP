classdef compute
    properties(Constant)
        subjids = {'MBC','MG','RC','WYZ','XLM','YC','YL','YMH','YZ','GG','SJ','XZ','TQ','AK'};
        model_names = {'opt','sum','max','min','corr','max2','min_max','min_corr','max_corr','min_max_corr',...
            'sum_erf','sum_min','sum_max','sum_corr','sum_min_max','sum_min_corr','sum_max_corr',...
            'sum_x_min', 'sum_x_max','sum_x_corr','sum_x_min_max','sum_x_min_corr','sum_x_max_corr','sum_x_min_max_corr', 'sign'...
            'random','err','sample','flex2','sump'};
        invalid_idx=[2,3,4,5,6,25,26,27,28,29];
        % shared directory
        dirs = 'real_data_results/fix_s/';
        dirs2 = 'raw_tables/';
        dirs3 = 'real_data_results/flex_s/'
        dirs_fake1 = 'fake_data_results/fix_s/';
        dirs_fake2 = 'fake_data_results/flex_s/';
        dirs_fakedata1 = 'fake_data/fix_s/';
        dirs_fakedata2 = 'fake_data/flex_s';
        
        % range for fake data parameters
        lambda_vec = linspace(0.001,0.3,31);
        guess_vec = linspace(0,1,51);
        lambda_fake = [0.01,0.05];
        guess_fake = [0,0.2];
        bin = linspace(-13,13,10);
        s_std = 9.0593;
        setsize = 3;
        sT_vec = linspace(-20,20,100);
        sD_vec = linspace(-20,20,100);
    end
    methods(Static)
        function obs_response = decision_rule(model_idx,x,lambda,nTrials,std_s)
            sigma = 1/sqrt(lambda);
            if ~exist('std_s','var')
                std_s = compute.s_std;
            end
            lambda_s = 1/std_s^2;
            switch model_idx
                case 1 % optimal model
%                     x_c = x*lambda/(lambda + lambda_s);
%                     std_c = 1/sqrt(lambda + lambda_s);
                    term = zeros(size(x));
                    for jj = 1:compute.setsize
                        term(jj,:) = ((sum(x)-x(jj,:))*lambda).^2/((compute.setsize-1)*lambda + lambda_s)/2 - (sum(x.^2) - x(jj,:).^2)*lambda/2;
                    end
                    obs_response = sum(erf(x*lambda/sqrt(2*(lambda+lambda_s))).*exp(-x.^2/(sigma^2+std_s^2)/2).*exp(term));
                case 2 % sum model
                    obs_response = sum(x);

                case 3 % max model
                    [~,idx] = max(abs(x));
                    idx = sub2ind(size(x), idx, 1:nTrials);
                    obs_response = x(idx);
                case 4 % min model
                    [~,idx] = min(abs(x));
                    idx = sub2ind(size(x), idx, 1:nTrials);
                    obs_response = x(idx);
                case 5 % corr model
                    corr_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp = x(idx_temp~=jj,:);
                        corr_x(jj,:) = std(temp);        
                    end
                    [~,idx] = min(corr_x);
                    idx = sub2ind(size(x), idx, 1:nTrials);	
                    obs_response = x(idx);

                case 6 % max2 model
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp = x(idx_temp~=jj,:);
                        term_x(jj,:) = mean(temp);
                    end
                    [~,idx] = min(abs(term_x));
                    idx = sub2ind(size(x),idx,1:nTrials);
                    obs_response = x(idx);

                case 7 % min-max model
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp1 = x(idx_temp==jj,:);
                        temp2 = x(idx_temp~=jj,:);
                        term_x(jj,:) = temp1.^2/(sigma^2 + std_s^2) + mean(temp2).^2/(sigma^2/(compute.setsize-1) + std_s^2);
                    end
                    [~,idx] = min(term_x);
                    idx = sub2ind(size(x), idx, 1:nTrials);	
                    obs_response = x(idx);
                case 8 % min-corr model
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp1 = x(idx_temp==jj,:);
                        temp2 = x(idx_temp~=jj,:);
                        term_x(jj,:) = temp1.^2/(sigma^2 + std_s^2) + ...
                            var(temp2,1)*(compute.setsize-1)*lambda;
                    end
                    [~,idx] = min(term_x);
                    idx = sub2ind(size(x), idx, 1:nTrials);	
                    obs_response = x(idx);    
                case 9 % max-corr model
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp = x(idx_temp~=jj,:);
                        term_x(jj,:) = mean(temp).^2/(sigma^2/(compute.setsize-1) + std_s^2) + var(temp,1)*(compute.setsize-1)*lambda;
                    end
                    [~,idx] = min(term_x);
                    idx = sub2ind(size(x), idx, 1:nTrials);	
                    obs_response = x(idx);
                case 10 % min-max-corr model                
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp1 = x(idx_temp==jj,:);
                        temp2 = x(idx_temp~=jj,:);
                        term_x(jj,:) = temp1.^2/(sigma^2 + std_s^2) + mean(temp2).^2/(sigma^2/(compute.setsize-1) + std_s^2) + var(temp2,1)*(compute.setsize-1)*lambda;
                    end
                    [~,idx] = min(term_x);
                    idx = sub2ind(size(x), idx, 1:nTrials);	
                    obs_response = x(idx); 

                case 11 % sum of erf model
                    x_c = x*lambda/(lambda + lambda_s);
                    std_c = 1/sqrt(lambda + lambda_s);

                    obs_response = sum(1-2*normcdf(0,x_c, std_c));
                case 12 % sum of min model
                    x_c = x*lambda/(lambda + lambda_s);
                    std_c = 1/sqrt(lambda + lambda_s);
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp1 = x(idx_temp==jj,:);
                        term_x(jj,:) = temp1.^2/(sigma^2 + std_s^2);
                    end
                    obs_response = sum(exp(-term_x/2).*(1-2*normcdf(0,x_c,std_c)));
                case 13 % sum of max model
                    x_c = x*lambda/(lambda + lambda_s);
                    std_c = 1/sqrt(lambda + lambda_s);
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp2 = x(idx_temp~=jj,:);
                        term_x(jj,:) = mean(temp2).^2/(sigma^2/(compute.setsize-1) + std_s^2);
                    end
                    obs_response = sum(exp(-term_x/2).*(1-2*normcdf(0,x_c,std_c)));
                case 14 % sum of variation model
                    x_c = x*lambda/(lambda + lambda_s);
                    std_c = 1/sqrt(lambda + lambda_s);
                    corr_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp = x(idx_temp~=jj,:);
                        corr_x(jj,:) = var(temp,1)*(compute.setsize-1)*lambda;
                    end
                    obs_response = sum(exp(-corr_x/2).*(1-2*normcdf(0,x_c,std_c)));
                case 15 % sum of min max model
                    x_c = x*lambda/(lambda + lambda_s);
                    std_c = 1/sqrt(lambda + lambda_s);
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp1 = x(idx_temp==jj,:);
                        temp2 = x(idx_temp~=jj,:);
                        term_x(jj,:) = temp1.^2/(sigma^2 + std_s^2) + mean(temp2).^2/(sigma^2/(compute.setsize-1) + std_s^2);
                    end
                    obs_response = sum(exp(-term_x/2).*(1-2*normcdf(0,x_c,std_c)));
                case 16 % sum of min corr model
                    x_c = x*lambda/(lambda + lambda_s);
                    std_c = 1/sqrt(lambda + lambda_s);
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp1 = x(idx_temp==jj,:);
                        temp2 = x(idx_temp~=jj,:);
                        term_x(jj,:) = temp1.^2/(sigma^2 + std_s^2) + ...
                            var(temp2,1)*(compute.setsize-1)*lambda;
                    end
                    obs_response = sum(exp(-term_x/2).*(1-2*normcdf(0,x_c,std_c)));
                case 17 % sum of max corr model
                    x_c = x*lambda/(lambda + lambda_s);
                    std_c = 1/sqrt(lambda + lambda_s);
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp = x(idx_temp~=jj,:);
                        term_x(jj,:) = mean(temp).^2/(sigma^2/(compute.setsize-1) + std_s^2) + var(temp,1)*(compute.setsize-1)*lambda;
                    end
                    obs_response = sum(exp(-term_x/2).*(1-2*normcdf(0,x_c,std_c)));
                case 18 % sum_x_min model
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp1 = x(idx_temp==jj,:);
                        term_x(jj,:) = temp1.^2/(sigma^2 + std_s^2);
                    end
                    obs_response = sum(x.*exp(-term_x/2));
                case 19 % sum_x_max model
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp2 = x(idx_temp~=jj,:);
                        term_x(jj,:) = mean(temp2).^2/(sigma^2/(compute.setsize-1) + std_s^2);
                    end
                    obs_response = sum(x.*exp(-term_x/2));
                case 20 % sum_x_corr model
                    corr_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp = x(idx_temp~=jj,:);
                        corr_x(jj,:) = var(temp,1)*(compute.setsize-1)*lambda;
                    end
                    obs_response = sum(x.*exp(-corr_x/2));
                case 21 % sum_x_min_max model
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp1 = x(idx_temp==jj,:);
                        temp2 = x(idx_temp~=jj,:);
                        term_x(jj,:) = temp1.^2/(sigma^2 + std_s^2) + mean(temp2).^2/(sigma^2/(compute.setsize-1) + std_s^2);
                    end
                    obs_response = sum(x.*exp(-term_x/2));
                case 22 % sum_x_min_corr model
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp1 = x(idx_temp==jj,:);
                        temp2 = x(idx_temp~=jj,:);
                        term_x(jj,:) = temp1.^2/(sigma^2 + std_s^2) + ...
                            var(temp2,1)*(compute.setsize-1)*lambda;
                    end
                    obs_response = sum(x.*exp(-term_x/2));
                case 23 % sum_x_max_corr model
                    term_x = zeros(size(x));
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp = x(idx_temp~=jj,:);
                        term_x(jj,:) = mean(temp).^2/(sigma^2/(compute.setsize-1) + std_s^2) + var(temp,1)*(compute.setsize-1)*lambda;
                    end
                    obs_response = sum(x.*exp(-term_x/2));
                case 24 % sum_x_min_max_corr model
                    term = zeros(size(x)); % compute weight
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp1 = x(idx_temp==jj,:);
                        temp2 = x(idx_temp~=jj,:);
                        term(jj,:) = temp1.^2/(sigma^2 + std_s^2) + mean(temp2).^2/(sigma^2/(compute.setsize-1) + std_s^2) + var(temp2,1)*(compute.setsize-1)*lambda;
                    end
                    obs_response = sum(x.*exp(-term/2));
                case 25 % sign model, observer reports the orientation that has the least number of items
                    temp = sign(x);
                    sum_x = sum(temp);
                    obs_response = zeros(1,length(x));
                    % all the measurements have the same sign
                    obs_response(sum_x==compute.setsize) = 1;
                    obs_response(sum_x==-compute.setsize) = -1;
                    % reports the orientation that has the least number of
                    % items, this only works for 4 items.
                    obs_response(sum_x==2) = -1;
                    obs_response(sum_x==-2) = 1;                 
                case 26 % random model
                    obs_response = 0;
                case 27 % err model
                    term = zeros(size(x)); % compute weight
                    idx_temp = 1:compute.setsize;
                    for jj = 1:compute.setsize
                        temp1 = x(idx_temp==jj,:);
                        temp2 = x(idx_temp~=jj,:);
                        term(jj,:) = temp1.^2/(sigma^2 + std_s^2) + mean(temp2).^2/(sigma^2/(compute.setsize-1) + std_s^2) + var(temp2)*(compute.setsize-1)*lambda;
                    end
                    obs_response = sum(sign(x).*exp(-term));
                    
                case 28 % sampling model
                    x_c = x*lambda/(lambda + lambda_s);
                    std_c = 1/sqrt(lambda + lambda_s);
                    term = zeros(size(x));
                    for jj = 1:compute.setsize
                        term(jj,:) = ((sum(x)-x(jj,:))*lambda).^2/((compute.setsize-1)*lambda + lambda_s)/2 - (sum(x.^2) - x(jj,:).^2)*lambda/2;
                    end
                    term1 = sum((1-normcdf(0,x_c, std_c)).*exp(-x.^2/(sigma^2+std_s^2)/2).*exp(term));
                    term2 = sum((normcdf(0,x_c, std_c)).*exp(-x.^2/(sigma^2+std_s^2)/2).*exp(term));
                    term = term1./term2;
                    obs_response = 1-1./(term+1); % probability of reporting right on this trial
            end
        end
        function generate_prediction_table_cluster(subj_idx, lambda_idx, model_idx)
            for j = model_idx
                dirname = compute.model_names{j};
                disp(['Computing prediction table for model ' dirname '...'])
                tic
                for i = subj_idx
                    
                    %initialize the random number generator
                    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
                    sigma = 0;
                    subjid = compute.subjids{i};
                    
                    %% read data
                    [stimuli,response,~] = utils.readdata(subjid);
                    %% parameters settings
                    lambda = linspace(0.001,0.3,31);
                    guess  = linspace(0,1,51);
                    sigma = 1./sqrt(lambda);
                    trial_num_sim = 1000; % trial number of simulation
                    stimuli = [stimuli(:,1) repmat(stimuli(:,2),1,compute.setsize-1)];

                    %% calculate the probability of certain response on each trial and combination of pars
                    log_CPG_prediction = zeros(length(lambda), length(guess));

                    %% calculate prediction
                    for kk = lambda_idx
                        % to save the prediction of a certain combination of pars
                        CP_prediction_temp = zeros(1, length(stimuli));
                        if ~exist([compute.dirs compute.dirs2 dirname],'dir')
                            mkdir([compute.dirs compute.dirs2 dirname]);
                        end
                        fname = [compute.dirs compute.dirs2 dirname '/' subjid '_' num2str(kk) '.mat'];
%                         if exist(fname,'file');
%                             continue
%                         end
                        for ii = 1:length(stimuli)
                           % generate noise samples
                            noiseMat = normrnd(zeros(compute.setsize,trial_num_sim), sigma(kk));
                            x = repmat(stimuli(ii,:),trial_num_sim,1);
                            x = x' + noiseMat;
                            
                            % decision rule
                            obs_response = compute.decision_rule(j,x,lambda(kk),trial_num_sim);
                            % calcualte the probability of reporting right
                            if j==28
                                CP_prediction_temp(ii) = mean(obs_response);
                            elseif j==26
                                CP_prediction_temp(ii) = 0.5;
                            else
                                CP_prediction_temp(ii) = (sum(obs_response>0) + .5*sum(obs_response==0))/trial_num_sim;
                            end
                        end
                        % generate the table for the CPG model
                        CP_prediction_temp_expanded = repmat(CP_prediction_temp, length(guess),1);
                        temp = bsxfun(@times, (1-guess'), CP_prediction_temp_expanded);
                        CPG_prediction_temp = bsxfun(@plus, temp, .5*guess');

                        CPG_prediction_temp(:,response==-1) = 1 - CPG_prediction_temp(:,response==-1);
                        CPG_prediction_temp(CPG_prediction_temp==1) = 1 - 1/trial_num_sim/10;
                        CPG_prediction_temp(CPG_prediction_temp==0) = 1/trial_num_sim/10;
                        % take log and sum over trials;
                        log_CPG_prediction(kk,:) = squeeze(sum(log(CPG_prediction_temp),2));
                        save(fname,'lambda','log_CPG_prediction','lambda_idx');
                        log_CPG_prediction = zeros(length(lambda),length(guess));
                        
                    end
                    
                    if ~exist([compute.dirs dirname '/pars'], 'dir');
                        mkdir([compute.dirs dirname '/pars']);
                    end
                    save([compute.dirs dirname '/pars/' subjid '_pars'], 'lambda', 'guess');
                end
            toc    
            end
            
        end
        function generate_prediction_table_full(subj_idx,model_idx)
            % create a full table of prediction from model 
            for j = model_idx
                dirname = compute.model_names{j};
                disp(['Computing prediction table for model ' dirname '...'])
                tic
                for i = subj_idx
                    
                    %initialize the random number generator
                    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
                    sigma = 0;
                    subjid = compute.subjids{i};
                    
                    %% read data
                    [stimuli,response,~] = utils.readdata(subjid);
                    %% parameters settings
                    lambda = linspace(0.001,0.3,31);
                    guess  = linspace(0,1,51);
                    sigma = 1./sqrt(lambda);
                    trial_num_sim = 1000; % trial number of simulation
                    stimuli = [stimuli(:,1) repmat(stimuli(:,2),1,compute.setsize-1)];

                    %% calculate the probability of certain response on each trial and combination of pars
                    predictionMat = zeros(length(lambda), length(guess),length(stimuli));

                    %% calculate prediction
                    for kk = 1:length(lambda)
                        % to save the prediction of a certain combination of pars
                        CP_prediction_temp = zeros(1, length(stimuli));
                        if ~exist([compute.dirs compute.dirs2 dirname],'dir')
                            mkdir([compute.dirs compute.dirs2 dirname]);
                        end
    
                        for ii = 1:length(stimuli)
                           % generate noise samples
                            noiseMat = normrnd(zeros(compute.setsize,trial_num_sim), sigma(kk));
                            x = repmat(stimuli(ii,:),trial_num_sim,1);
                            x = x' + noiseMat;
                            
                            % decision rule
                            obs_response = compute.decision_rule(j,x,lambda(kk),trial_num_sim);
                            % calcualte the probability of reporting right
                            if j==28
                                CP_prediction_temp(ii) = mean(obs_response);
                            elseif j==26
                                CP_prediction_temp(ii) = 0.5;
                            else
                                CP_prediction_temp(ii) = (sum(obs_response>0) + .5*sum(obs_response==0))/trial_num_sim;
                            end
                        end
                        % generate the table for the CPG model
                        CP_prediction_temp_expanded = repmat(CP_prediction_temp, length(guess),1);
                        temp = bsxfun(@times, (1-guess'), CP_prediction_temp_expanded);
                        CPG_prediction_temp = bsxfun(@plus, temp, .5*guess');

                        CPG_prediction_temp(:,response==-1) = 1 - CPG_prediction_temp(:,response==-1);
                        CPG_prediction_temp(CPG_prediction_temp==1) = 1 - 1/trial_num_sim/10;
                        CPG_prediction_temp(CPG_prediction_temp==0) = 1/trial_num_sim/10;

                        predictionMat(kk,:,:) = CPG_prediction_temp;
                    end
                    % compute likelihood of a model for each trial
                    predictionMat = reshape(predictionMat,[length(lambda)*length(guess),length(stimuli)]);
                    logMat = log(predictionMat);
                    L_max = max(predictionMat);
                    L = bsxfun(@minus,logMat,L_max);
                    prediction = L_max + log(sum(exp(L))) - log(length(lambda)*length(guess));
                    prediction = exp(prediction);
                    
                    if ~exist([compute.dirs dirname '/tables_full'], 'dir');
                        mkdir([compute.dirs dirname '/tables_full']);
                    end
                    save([compute.dirs dirname '/tables_full/' subjid '.mat'], 'prediction');
                end
            toc    
            end
            
        end
        function generate_prediction_table_cluster_flex_s(subj_idx, lambda_idx, model_idx)
            for i = subj_idx
                for j = model_idx
                tic
                    if ~isempty(find(compute.invalid_idx==j))
                        disp('This model is invalid for this computation');
                        continue
                    end
                    %initialize the random number generator
                    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
                    sigma = 0;
                    subjid = compute.subjids{i};
                    dirname = compute.model_names{j};
                    %% read data
                    [stimuli,response,~] = utils.readdata(subjid);
                    %% parameters settings
                    lambda = linspace(0.001,0.3,31);
                    std_s_vec = linspace(1,25,31);
                    guess  = linspace(0,1,51);
                    sigma = 1./sqrt(lambda);
                    trial_num_sim = 1000; % trial number of simulation
                    stimuli = [stimuli(:,1) repmat(stimuli(:,2),1,compute.setsize-1)];

                    %% calculate the probability of certain response on each trial and combination of pars
                    log_CPG_prediction = zeros(length(lambda), length(std_s_vec), length(guess));

                    %% calculate prediction
                    for kk = lambda_idx
                        for jj = 1:length(std_s_vec)
                            % to save the prediction of a certain combination of pars
                            CP_prediction_temp = zeros(1, length(stimuli));
                            if ~exist([compute.dirs3 compute.dirs2 dirname],'dir')
                                mkdir([compute.dirs3 compute.dirs2 dirname]);
                            end
                            fname = [compute.dirs3 compute.dirs2 dirname '/' subjid '_' num2str(kk) '.mat'];
%                             if exist(fname,'file');
%                                 continue
%                             end
                            for ii = 1:length(stimuli)
                               % generate noise samples
                                noiseMat = normrnd(zeros(compute.setsize,trial_num_sim), sigma(kk));
                                x = repmat(stimuli(ii,:),trial_num_sim,1);
                                x = x' + noiseMat;

                                % decision rule
                                obs_response = compute.decision_rule(j,x,lambda(kk),trial_num_sim,std_s_vec(jj));
                                % calcualte the probability of reporting right

                                CP_prediction_temp(ii) = (sum(obs_response>0) + .5*sum(obs_response==0))/trial_num_sim;

                            end
                        
                            % generate the table for the CPG model
                            CP_prediction_temp_expanded = repmat(CP_prediction_temp, length(guess),1);
                            temp = bsxfun(@times, (1-guess'), CP_prediction_temp_expanded);
                            CPG_prediction_temp = bsxfun(@plus, temp, .5*guess');

                            CPG_prediction_temp(:,response==-1) = 1 - CPG_prediction_temp(:,response==-1);
                            CPG_prediction_temp(CPG_prediction_temp==1) = 1 - 1/trial_num_sim/10;
                            CPG_prediction_temp(CPG_prediction_temp==0) = 1/trial_num_sim/10;
                            % take log and sum over trials;
                            log_CPG_prediction(kk,jj,:) = squeeze(sum(log(CPG_prediction_temp),2));
                              
                        end
                        save(fname,'lambda','log_CPG_prediction','lambda_idx');
                            log_CPG_prediction = zeros(length(lambda),length(std_s_vec),length(guess));
                    end
                    if ~exist([compute.dirs3 dirname '/pars'], 'dir');
                        mkdir([compute.dirs3 dirname '/pars']);
                    end
                    toc
                    save([compute.dirs3 dirname '/pars/' subjid '_pars'], 'lambda', 'std_s_vec', 'guess');
                end
            end
            
        end
        function generate_prediction_table_cluster_x_p(subj_idx, lambda_idx)
            for i = subj_idx
                tic
                    %initialize the random number generator
                    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
                    sigma = 0;
                    subjid = compute.subjids{i};
                    dirname = 'sump';
                    %% read data
                    [stimuli,response,~] = utils.readdata(subjid);
                    %% parameters settings
                    lambda = linspace(0.001,0.3,31);
                    p = linspace(0,5,31);
                    guess  = linspace(0,1,51);
                    sigma = 1./sqrt(lambda);
                    trial_num_sim = 1000; % trial number of simulation
                    stimuli = [stimuli(:,1) repmat(stimuli(:,2),1,compute.setsize-1)];

                    %% calculate the probability of certain response on each trial and combination of pars
                    log_CPG_prediction = zeros(length(lambda), length(p), length(guess));

                    %% calculate prediction
                    for kk = lambda_idx
                        for jj = 1:length(p)
                            % to save the prediction of a certain combination of pars
                            CP_prediction_temp = zeros(1, length(stimuli));
                            if ~exist([compute.dirs compute.dirs2 dirname],'dir')
                                mkdir([compute.dirs compute.dirs2 dirname]);
                            end
                            fname = [compute.dirs compute.dirs2 dirname '/' subjid '_' num2str(kk) '.mat'];
%                             if exist(fname,'file');
%                                 continue
%                             end
                            for ii = 1:length(stimuli)
                               % generate noise samples
                                noiseMat = normrnd(zeros(compute.setsize,trial_num_sim), sigma(kk));
                                x = repmat(stimuli(ii,:),trial_num_sim,1);
                                x = x' + noiseMat;

                                % decision rule
                                obs_response = sum(x.^p(jj));
                                % calcualte the probability of reporting right

                                CP_prediction_temp(ii) = (sum(obs_response>0) + .5*sum(obs_response==0))/trial_num_sim;

                            end
                        
                            % generate the table for the CPG model
                            CP_prediction_temp_expanded = repmat(CP_prediction_temp, length(guess),1);
                            temp = bsxfun(@times, (1-guess'), CP_prediction_temp_expanded);
                            CPG_prediction_temp = bsxfun(@plus, temp, .5*guess');

                            CPG_prediction_temp(:,response==-1) = 1 - CPG_prediction_temp(:,response==-1);
                            CPG_prediction_temp(CPG_prediction_temp==1) = 1 - 1/trial_num_sim/10;
                            CPG_prediction_temp(CPG_prediction_temp==0) = 1/trial_num_sim/10;
                            % take log and sum over trials;
                            log_CPG_prediction(kk,jj,:) = squeeze(sum(log(CPG_prediction_temp),2));
                              
                        end
                        save(fname,'lambda','log_CPG_prediction','lambda_idx');
                            log_CPG_prediction = zeros(length(lambda),length(p),length(guess));
                    end
                    if ~exist([compute.dirs dirname '/pars'], 'dir');
                        mkdir([compute.dirs dirname '/pars']);
                    end
                    toc
                    save([compute.dirs dirname '/pars/' subjid '_pars'], 'lambda', 'p', 'guess');
            end
            
        end
        function generate_prediction_table_cluster_fake_data(lambda_idx,model_idx)
                           
            for nn = model_idx
                for mm = lambda_idx
                    %initialize the random number generator
                    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
                    dirname = compute.model_names{nn};
                    %% read data
                    load('fake_data_stimuli.mat');
                    stimuli = stimuli_fake';
                    trial_num = length(stimuli);
                    
                    lambda = compute.lambda_vec;

                    trial_num_sim = 1000; % trial number of simulation
                    stimuli = [stimuli(:,1) repmat(stimuli(:,2),1,compute.setsize-1)];

                    %% calculate the probability of certain response on each trial and combination of pars
                    CP_prediction_r = zeros(trial_num, length(lambda));
                  %% CP model
                    fname = [compute.dirs_fake1 dirname '/raw_tables/prediction_table_' num2str(mm) '.mat'];
                    for kk = mm
                        for ii = 1:length(stimuli)
                            % generate noise samples
                            noiseMat = normrnd(zeros(compute.setsize,trial_num_sim), 1/sqrt(lambda(kk)));
                            x = repmat(stimuli(ii,:),trial_num_sim,1);
                            x = x' + noiseMat;
                            
                            % decision rule
                            obs_response = compute.decision_rule(nn,x,lambda(kk),trial_num_sim);
                                                 
                            % calculate the probability of reporting right
                            if nn==26
                                CP_prediction_r(ii,kk) = mean(obs_response);
                            else
                                CP_prediction_r(ii,kk) = (sum(obs_response>0) + .5*sum(obs_response==0))/trial_num_sim; 
                            end
                        end
                    end
                    if ~exist([compute.dirs_fake1 dirname '/raw_tables/'],'dir');
                        mkdir([compute.dirs_fake1 dirname '/raw_tables/']);
                    end
                    save(fname,'lambda','CP_prediction_r','lambda_idx');

                end
            end
        end
        function generate_prediction_table_grid(lambda_idx,model_idx)
            sT = compute.sT_vec;
            sD = compute.sD_vec;
            for nn = model_idx
                dirname = compute.model_names{nn};
                tic
                disp(['Computing prediction table for model ' dirname '...'])
                for mm = lambda_idx
                    %initialize the random number generator
                    RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
                    
                    %% read data
                    
                    lambda = compute.lambda_vec;

                    trial_num_sim = 1000; % trial number of simulation

                    %% calculate the probability of certain response on each trial and combination of pars
                    CP_prediction_r = zeros(length(sT),length(sD), length(lambda));
                  %% CP model
                    fname = [compute.dirs dirname '/tables_grid/prediction_table_' num2str(mm) '.mat'];
                    
                    for ii = 1:length(sT)
                        for jj = 1:length(sD)
                            % generate noise samples
                            stimuli = [sT(ii) repmat(sD(jj),1,compute.setsize-1)];
                            noiseMat = normrnd(zeros(compute.setsize,trial_num_sim), 1/sqrt(lambda(mm)));
                            x = repmat(stimuli,trial_num_sim,1);
                            x = x' + noiseMat;

                            % decision rule
                            obs_response = compute.decision_rule(nn,x,lambda(mm),trial_num_sim);

                            % calculate the probability of reporting right
                            if nn==28
                                CP_prediction_r(ii,jj,mm) = mean(obs_response);
                            else
                                CP_prediction_r(ii,jj,mm) = (sum(obs_response>0) + .5*sum(obs_response==0))/trial_num_sim; 
                            end
                        end
                    end
                   
                    if ~exist([compute.dirs dirname '/tables_grid/'],'dir');
                        mkdir([compute.dirs dirname '/tables_grid/']);
                    end
                    save(fname,'lambda','CP_prediction_r','lambda_idx');

                end
                toc
            end
        end
        function combine_tables(subj_idx, model_idx, flex_s)
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            for ii = subj_idx
                for jj = model_idx
                    subjid = compute.subjids{ii};
                    dirname = compute.model_names{jj};
                    if flex_s
                        dirname2 = [compute.dirs3 compute.dirs2 dirname];
                        dirname3 = [compute.dirs3 dirname];
                    else
                        dirname2 = [compute.dirs compute.dirs2 dirname];
                        dirname3 = [compute.dirs dirname];
                    end
                    load([dirname2 '/' subjid '_1.mat']);
                    log_CPG_prediction_all = zeros(size(log_CPG_prediction));
                    for kk = 1:length(lambda)
                        load([dirname2 '/' subjid '_' num2str(kk) '.mat']);
                        log_CPG_prediction_all = log_CPG_prediction_all + log_CPG_prediction;            
                    end
                    log_CPG_prediction = log_CPG_prediction_all;
                    if ~exist([dirname3 '/tables'], 'dir')
                        mkdir([dirname3 '/tables']);
                    end
                    save([dirname3 '/tables/' compute.subjids{ii} '.mat'],'log_CPG_prediction');
                end
            end
        end
        function combine_tables_fake(model_idx,flex_s)
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            if flex_s
                common_dir = compute.dirs_fake2;
            else
                common_dir = compute.dirs_fake1;
            end
            for kk = model_idx
                        
                dirname = compute.model_names{kk};
                load([common_dir dirname '/raw_tables/prediction_table_1.mat']);
                CP_prediction_r_all = zeros(size(CP_prediction_r));

                for mm=1:length(lambda)
                    load([common_dir dirname '/raw_tables/prediction_table_' num2str(mm) '.mat']);
                    CP_prediction_r_all = CP_prediction_r_all + CP_prediction_r;            
                end

                CP_prediction_r = CP_prediction_r_all;
                save([common_dir dirname '/tables.mat'], 'CP_prediction_r');
           end
                
        end
        function combine_tables_grid(model_idx)
            
            for kk = model_idx
                        
                dirname = compute.model_names{kk};
                load([compute.dirs dirname '/tables_grid/prediction_table_1.mat']);
                CP_prediction_r_all = zeros(size(CP_prediction_r));

                for mm=1:length(lambda)
                    load([compute.dirs dirname '/tables_grid/prediction_table_' num2str(mm) '.mat']);
                    CP_prediction_r_all = CP_prediction_r_all + CP_prediction_r;            
                end

                CP_prediction_r = CP_prediction_r_all;
                save([compute.dirs dirname '/tables_grid.mat'], 'CP_prediction_r');
           end
                
        end
        function p_right_bin(subj_idx,half)
            % half secifies which part of data to use
            % 0: whole data set
            % 1: data with odd indices
            % 2: data with even indices
            % modfied 2015-09-14 SS
            for kk = subj_idx
                subjid = compute.subjids{kk};
                [stimuli,response] = utils.readdata(subjid,half);
                
                load('bins_new');
                [cnt,cnt_r,p_right] = utils.compute_cnts(stimuli,response,bins);
               
                if ~exist([compute.dirs 'data_bin_' num2str(half)],'dir')
                    mkdir([compute.dirs 'data_bin_' num2str(half)])
                end
                save([compute.dirs 'data_bin_' num2str(half) '/' subjid '.mat'],'p_right','cnt_r','cnt')
            end
        end
        function fit_pars(subj_idx,model_idx,varargin)
            % compute fit parameters from table
            % modified to work for both real data and fake data
            % 15-07-08 SS
            [flex_s, ~, ~, common_dir,subj_ids] = utils.parse_inputs(subj_idx, varargin);
            
            for ii = 1:length(subj_ids)
                for jj = model_idx
                    subjid = subj_ids{ii};
                    dirname = compute.model_names{jj};
                    load([common_dir dirname '/tables/' subjid '.mat']);
                    [~,I] = max(log_CPG_prediction(:));
                    lambda = compute.lambda_vec;
                    guess = compute.guess_vec;
                    if flex_s
                        [idx1,idx2,idx3] = ind2sub(size(log_CPG_prediction),I);
                        CPG_lambda_hat = lambda(idx1);
                        CPG_std_s_hat = std_s_vec(idx2); 
                        CPG_guess_hat = guess(idx3);
                    else
                        [idx1,idx2] = ind2sub(size(log_CPG_prediction),I);
                        CPG_lambda_hat = lambda(idx1);
                        CPG_guess_hat = guess(idx2);
                    end
                    if ~exist([common_dir dirname '/fit_pars'], 'dir')
                        mkdir([common_dir dirname '/fit_pars']);
                    end
                    if flex_s
                        save([common_dir dirname '/fit_pars/' subjid],'CPG_lambda_hat','CPG_guess_hat','CPG_std_s_hat');
                    else
                        save([common_dir dirname '/fit_pars/' subjid],'CPG_lambda_hat','CPG_guess_hat');
                    end
                end
            end
        end     
        function evi(subj_idx, model_idx,varargin)
            % compute bmc, bic, aic, aic, aicc of a certain model for a
            % certain subject
            % For real data, subj_idx is the idx of subjects in
            % compute.subjids. For fake data, subj_idx is the model idx to
            % generate fake data
            % model_idx is the idx to fit the data (real or fake)
            % modified to fit both real data and fake data
            % 15-07-08 SS
            
            
            % parse the conditions: fake or real data, whether has flex s
            [flex_s, npars, nTrials, common_dir, subj_ids] = utils.parse_inputs(subj_idx, varargin);
            
            
            for ii = 1:length(subj_ids)
                for jj = model_idx
                    subjid = subj_ids{ii};
                    dirname = compute.model_names{jj};
                    if ~isempty(find(compute.invalid_idx==jj)) && flex_s
                        disp('This model is invalid for this computation');
                        continue
                    end
                    load([common_dir dirname '/tables/' subjid '.mat']);
                    % bmc
                    L_max_CPG = max(log_CPG_prediction(:));                    
                    L_CPG = log_CPG_prediction - L_max_CPG;
                    bmc = L_max_CPG + log(sum(exp(L_CPG(:)))) - log(numel(L_CPG));
                    % bic
                    bic = L_max_CPG - 0.5*npars*log(nTrials);
                    % aic
                    aic = L_max_CPG - npars;
                    % aicc
                    aicc = aic - npars*(npars+1)/(nTrials-npars-1);
                    if ~exist([common_dir dirname '/evi'], 'dir')
                        mkdir([common_dir dirname '/evi']);
                    end
                    save([common_dir dirname '/evi/' subjid], 'aic', 'bic', 'bmc','aicc');
                end
            end
        end 
        function evi_bin(subj_idx, model_idx)
            % compute LML of each bin
            for ii = 1:length(subj_idx)
                subjid = compute.subjids{subj_idx(ii)};
                stimuli = utils.readdata(subjid);
                target_stimuli = stimuli(:,1);
                dist_stimuli = stimuli(:,2);
                load('bins_new');
                idx1 = interp1(bins,1:length(bins),target_stimuli, 'nearest','extrap');
                idx2 = interp1(bins,1:length(bins),dist_stimuli, 'nearest','extrap');
                pred = zeros(length(bins),length(bins));
                cnt = zeros(size(pred));
                for mm = 1:length(model_idx)
                    dirname = compute.model_names{model_idx(mm)};
                    load([compute.dirs dirname '/tables_full/' subjid '.mat'])
                    for jj = 1:length(bins)
                        for kk = 1:length(bins)
                            pred_rel = prediction(idx1==jj & idx2==kk);
                            pred(jj,kk) = mean(pred_rel);
                            cnt(jj,kk) = length(pred_rel);
                        end
                    end
                    log_pred = log(pred).*cnt;
                    bmc = sum(log_pred(:));
                    if ~exist([compute.dirs dirname '/evi_bin'],'dir')
                        mkdir([compute.dirs dirname '/evi_bin'])
                    end
                    save([compute.dirs dirname '/evi_bin/' subjid], 'bmc','pred');
                end
            end
        end
        function fit_predictions(subj_idx, model_idx, varargin)
            % parse the conditions: fake or real data, whether has flex s
            [flex_s, ~, ~, common_dir, subj_ids] = utils.parse_inputs(subj_idx, varargin);
            
            for mm = 1:length(subj_ids)
                for nn = model_idx
                    % calculate fit predictions for models and p_right for real data
                    dirname = compute.model_names{nn};
                    subjid = subj_ids{mm};
                    %% load data
                    load([common_dir dirname '/fit_pars/' subjid]);
                    
                    if strcmp(common_dir, 'real_data_results/fix_s/')
                        [stimuli,~,~] = utils.readdata(subjid);
                    else
                        load('fake_data_stimuli','stimuli_fake');
                        stimuli = stimuli_fake';
                    end
                    trial_num = length(stimuli);
                    trial_num_sim = 1000;
                    %% pars settings
                    stimuli = [stimuli(:,1) repmat(stimuli(:,2),1,compute.setsize-1)];             
                    %% CP & CPG models
                    % loop over real trials, calculate the prediction for each stimulus
                    for i = 1:trial_num
                        noiseMat_CPG = normrnd(0,repmat(1/sqrt(CPG_lambda_hat),compute.setsize,trial_num_sim));
                        
                        %% CPG
                        % compute predictions
                        x = repmat(stimuli(i,:),trial_num_sim,1);
                        x = x' + noiseMat_CPG;
                        % decision rule
                        if flex_s
                            obs_response = compute.decision_rule(nn,x,CPG_lambda_hat,trial_num_sim,CPG_std_s_hat);
                        else
                            obs_response = compute.decision_rule(nn,x,CPG_lambda_hat,trial_num_sim);
                        end
                        % the prediction of the probability of reporting right
                        if nn==26
                            CPG_prediction(i) = mean(obs_response);
                        else
                            CPG_prediction(i) = (sum(obs_response>0) + .5*sum(obs_response==0))/trial_num_sim;
                        end

                    end
                    CPG_prediction = CPG_prediction*(1-CPG_guess_hat) + .5*CPG_guess_hat;

                    %% save result
                    if ~exist([common_dir dirname '/fit_results'], 'dir')
                        mkdir([common_dir dirname '/fit_results']);
                    end
                    save([common_dir dirname '/fit_results/' subjid], 'CPG_prediction');
                end
            end
        end   
        function fit_predictions_1d(subj_idx, model_idx,flex_s)
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            if flex_s
                common_dir = compute.dirs3;
            else
                common_dir = compute.dirs;
            end
            for mm = subj_idx
                for nn = model_idx
                    % calculate fit predictions for models and p_right for real data
                    dirname = compute.model_names{nn};
                    subjid = compute.subjids{mm};
                    %% load data
                    load([common_dir dirname '/fit_results/' subjid]);
                    [stimuli response] = utils.readdata(subjid);
                    target_stimuli = stimuli(:,1);
                    dist_stimuli = stimuli(:,2);
                    %% pars settings
                    bins_target = -15:3:15;
                    bins_dist = [-10,0,10];
                    
                    %% calculate the predictions of the fitted pars
                    prediction = zeros(length(bins_target), length(bins_dist));
            
                    %% CP & CPG models                  
                    idx1 = interp1(bins_target,1:length(bins_target),target_stimuli, 'nearest','extrap');
                    idx2 = interp1(bins_dist,1:length(bins_dist),dist_stimuli, 'nearest','extrap');
                    for ii = 1:length(bins_target)
                        for jj = 1:length(bins_dist)
                            CPG_prediction_rel = CPG_prediction(idx1==ii & idx2==jj);
                            prediction(ii,jj) = mean(CPG_prediction_rel);
                        end
                    end

                    %% calculate the real probability of reporting right
                    target_stimuli_right = target_stimuli(response==1);
                    dist_stimuli_right = dist_stimuli(response==1);

                    cnt_r = hist2d([target_stimuli_right,dist_stimuli_right], bins_target, bins_dist);
                    cnt = hist2d([target_stimuli,dist_stimuli],bins_target,bins_dist);
                    p_right = cnt_r./cnt;
                    cntVec = cnt;

                    %% save result                  
                    if ~exist([common_dir dirname '/fit_results_1d'], 'dir')
                        mkdir([common_dir dirname '/fit_results_1d']);
                    end
                    save([common_dir dirname '/fit_results_1d/' subjid], 'p_right','prediction','bins_target','bins_dist', 'cntVec');
                end
            end
        end   
        function fit_predictions_1d_dist(subj_idx, model_idx, flex_s)
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            if flex_s
                common_dir = compute.dirs3;
            else
                common_dir = compute.dirs;
            end
            for mm = subj_idx
                for nn = model_idx
                    % calculate fit predictions for models and p_right for real data
                    dirname = compute.model_names{nn};
                    subjid = compute.subjids{mm};
                    %% load data
                    load([common_dir dirname '/fit_results/' subjid]);
                    [stimuli response] = utils.readdata(subjid);
                    target_stimuli = stimuli(:,1);
                    dist_stimuli = stimuli(:,2);
                    %% pars settings
                    bins_dist = -15:3:15;
                    bins_target = [-10,0,10];

                    %% calculate the predictions of the fitted pars
                    prediction = zeros(length(bins_target), length(bins_dist));
                    
                    idx1 = interp1(bins_target,1:length(bins_target),target_stimuli, 'nearest','extrap');
                    idx2 = interp1(bins_dist,1:length(bins_dist),dist_stimuli, 'nearest','extrap');
                    for ii = 1:length(bins_target)
                        for jj = 1:length(bins_dist)
                            CPG_prediction_rel = CPG_prediction(idx1==ii & idx2==jj);
                            prediction(ii,jj) = mean(CPG_prediction_rel);
                        end
                    end

                    %% calculate the real probability of reporting right
                    target_stimuli_right = target_stimuli(response==1);
                    dist_stimuli_right = dist_stimuli(response==1);

                    cnt_r = hist2d([target_stimuli_right,dist_stimuli_right], bins_target, bins_dist);
                    cnt = hist2d([target_stimuli,dist_stimuli],bins_target,bins_dist);
                    p_right = cnt_r./cnt;
                    cntVec = cnt;

                    %% save results
                    if ~exist([common_dir dirname '/fit_results_1d_dist'], 'dir')
                        mkdir([common_dir dirname '/fit_results_1d_dist']);
                    end
                    save([common_dir dirname '/fit_results_1d_dist/' subjid], 'p_right','prediction','bins_target','bins_dist', 'cntVec');
                end
            end
        end  
        function fit_predictions_2d(subj_idx, model_idx, flex_s)
            if ~exist('flex_s', 'var')
                flex_s=0;
            end
            if flex_s
                common_dir = compute.dirs3;
            else
                common_dir = compute.dirs;
            end
            for mm = subj_idx
                for nn = model_idx
                    dirname = compute.model_names{nn};
                    subjid = compute.subjids{mm};
                    load([compute.dirs dirname '/fit_results/' subjid]);
                    [stimuli,response] = utils.readdata(subjid);
                    target_stimuli = stimuli(:,1);
                    dist_stimuli = stimuli(:,2);
                    %% pars settings
                    bins = compute.bin;
                    %% calculate the predictions of the fitted pars
                    prediction = zeros(length(bins), length(bins));
                    idx1 = interp1(bins,1:length(bins),target_stimuli, 'nearest','extrap');
                    idx2 = interp1(bins,1:length(bins),dist_stimuli, 'nearest','extrap');
                    
                    for ii = 1:length(bins)
                        for jj = 1:length(bins)
                            CPG_prediction_rel = CPG_prediction(idx1==ii & idx2==jj);
                            prediction(ii,jj) = mean(CPG_prediction_rel);
                        end
                    end

                    %% calculate the real probability of reporting right
                    target_stimuli_right = target_stimuli(response==1);
                    dist_stimuli_right = dist_stimuli(response==1);

                    cnt_r = hist2d([target_stimuli_right,dist_stimuli_right], bins, bins);
                    cnt = hist2d([target_stimuli,dist_stimuli],bins,bins);
                    p_right = cnt_r./cnt;
                    cntVec = cnt;

                    %% save results
                    
                    if ~exist([compute.dirs dirname '/fit_results_2d'], 'dir')
                        mkdir([compute.dirs dirname '/fit_results_2d']);
                    end
                    save([compute.dirs dirname '/fit_results_2d/' subjid], 'p_right','prediction','bins', 'cntVec');
                end
            end
        end
        function fit_predictions_2d2(subj_idx, model_idx,flex_s)
            % calculate the fit predictions in each grid of quantile.
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            
            if flex_s
                common_dir = compute.dirs3;
            else
                common_dir = compute.dirs;
            end
            % plot psychometric color map with new bins
            for mm = subj_idx
                for nn = model_idx
                    dirname = compute.model_names{nn};
                    subjid = compute.subjids{mm};
                    load([common_dir dirname '/fit_results/' subjid]);
                    [stimuli,response] = utils.readdata(subjid);
                    target_stimuli = stimuli(:,1);
                    dist_stimuli = stimuli(:,2);
                    
                    % load bins
                    load('bins_quantile');
                    
                    %% calculate the predictions of the fitted pars
                    prediction = zeros(length(bins), length(bins));
                    
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
                            CPG_prediction_rel = CPG_prediction(idx1==ii & idx2==jj);
                            prediction(ii,jj) = mean(CPG_prediction_rel);
                            cnt(ii,jj) = length(CPG_prediction_rel);
                            temp = target_stimuli(idx1_right==ii & idx2_right==jj);
                            cnt_r(ii,jj) = length(temp);
                        end
                    end
                    %% calculate the real probability of reporting right
                    p_right = cnt_r./cnt;
                    cntVec = cnt;

                    %% save results
                    
                    if ~exist([common_dir dirname '/fit_results_2d2'], 'dir')
                        mkdir([common_dir dirname '/fit_results_2d2']);
                    end
                    save([common_dir dirname '/fit_results_2d2/' subjid], 'p_right','prediction','bins', 'cntVec');
                end
            end
        end
        function fit_performance_2d2(subj_idx, model_idx)
            % plot psychometric color map with new bins
            for mm = subj_idx
                for nn = model_idx
                    dirname = compute.model_names{nn};
                    subjid = compute.subjids{mm};
                    load([compute.dirs dirname '/fit_results/' subjid]);
                    [stimuli,~,performance] = utils.readdata(subjid);
                    target_stimuli = stimuli(:,1);
                    dist_stimuli = stimuli(:,2);
                    
                    % load bins
                    load('bins_quantile');
                    

                    %% calculate the predictions of the fitted pars
                    prediction = zeros(length(bins), length(bins));
                    
                    % real response
                    target_stimuli_right = target_stimuli(performance==1);
                    dist_stimuli_right = dist_stimuli(performance==1);
                    idx1 = interp1(bins,1:length(bins),target_stimuli, 'nearest','extrap');
                    idx2 = interp1(bins,1:length(bins),dist_stimuli, 'nearest','extrap');
                    idx1_right = interp1(bins,1:length(bins),target_stimuli_right, 'nearest','extrap');
                    idx2_right = interp1(bins,1:length(bins),dist_stimuli_right, 'nearest','extrap');
                    cnt_r = zeros(length(bins), length(bins));
                    cnt = zeros(length(bins), length(bins));
                    for ii = 1:length(bins)
                        for jj = 1:length(bins)
                            CPG_prediction_rel = CPG_prediction(idx1==ii & idx2==jj);
                            prediction(ii,jj) = mean(CPG_prediction_rel);
                            cnt(ii,jj) = length(CPG_prediction_rel);
                            temp = target_stimuli(idx1_right==ii & idx2_right==jj);
                            cnt_r(ii,jj) = length(temp);
                        end
                    end

                    %% calculate the real probability of reporting right
                    p_right = cnt_r./cnt;
                    cntVec = cnt;

                    %% save results
                    if ~exist([compute.dirs dirname '/fit_results_2d2_perf'], 'dir')
                        mkdir([compute.dirs dirname '/fit_results_2d2_perf']);
                    end
                    save([compute.dirs dirname '/fit_results_2d2/perf' subjid], 'p_right','prediction','bins', 'cntVec');
                end
            end
        end
        function performance_curve(model_idx,sigma)
            for mm = model_idx
                % load stimuli
                tic
                load('stimuli.mat');
                stimuli = [sT; repmat(sD,compute.setsize-1,1)];
                performance = zeros(size(sigma));
                for kk = 1:length(sigma)
                    lambda = 1/sigma(kk)^2;
                    x = normrnd(stimuli, sigma(kk));
                     % decision rule
                    obs_response = compute.decision_rule(mm,x,lambda,length(stimuli));
                    obs_response(obs_response>0) = 1; obs_response(obs_response<=0) = -1;
                    
                    truth = sT;
                    truth(truth>0) = 1;
                    truth(truth<=0) = -1;
                    obs_response(obs_response~=truth) = 0;
                    obs_response(obs_response==truth) = 1;
                    performance(kk) = mean(obs_response);
                end
                dirname = compute.model_names{mm};
                if ~exist([compute.dirs dirname '/perf/'], 'dir');
                        mkdir([compute.dirs dirname '/perf/']);
                end
                save([compute.dirs dirname '/perf/performance.mat'  ], 'performance','sigma');
                toc
            end
        end
        function prediction_relative_to_opt(model_idx,sigma)
            load('stimuli.mat');
            std_s = compute.s_std;             lambda_s = 1/std_s^2;
            stimuli = [sT; repmat(sD,compute.setsize-1,1)];
            obs_response_opt = zeros(length(sigma),length(stimuli));
            for kk = 1:length(sigma)
                lambda = 1/sigma(kk)^2;
                x = normrnd(stimuli, sigma(kk));
                x_c = x*lambda/(lambda + lambda_s);
                std_c = 1/sqrt(lambda + lambda_s);
                term = zeros(size(x));
                
                for jj = 1:compute.setsize
                    term(jj,:) = ((sum(x)-x(jj,:))*lambda).^2/((compute.setsize-1)*lambda + lambda_s)/2 - (sum(x.^2) - x(jj,:).^2)*lambda/2;
                end
                obs_response = sum((1-2*normcdf(0,x_c, std_c)).*exp(-x.^2/(sigma(kk)^2+std_s^2)/2).*exp(term));
                obs_response(obs_response>0)=1; obs_response(obs_response<=0)=-1;
                obs_response_opt(kk,:) = obs_response;
            end
            for mm = model_idx               
                pred2opt = zeros(size(sigma));
                for kk = 1:length(sigma)
                    lambda = 1/sigma(kk)^2;
                    x = normrnd(stimuli, sigma(kk));
                     % decision rule
                    obs_response = compute.decision_rule(mm,x,lambda,length(stimuli));
                    obs_response(obs_response>0) = 1; obs_response(obs_response<=0) = -1;
                    
                    pred2opt(kk) = mean(obs_response==obs_response_opt(kk,:));
                end
                dirname = compute.model_names{mm};
                if ~exist([compute.dirs dirname '/pred2opt/'], 'dir');
                        mkdir([compute.dirs dirname '/pred2opt/']);
                end
                save([compute.dirs dirname '/pred2opt/pred2opt.mat'  ], 'pred2opt','sigma');
            end
        end
        function prediction_relative_to_opt2(model_idx)
            % given same measurements, compute the probability of giving
            % the same prediction comparing to optimal model
            
            sigma = 0.01:0.01:10;
            load('stimuli.mat');
            xMat = utils.generate_measurements(sigma);
            std_s = compute.s_std;             lambda_s = 1/std_s^2;
            
            obs_response_opt = zeros(length(sigma),length(xMat));
            for kk = 1:length(sigma)
                lambda = 1/sigma(kk)^2;
                x = squeeze(xMat(:,:,kk));
                x_c = x*lambda/(lambda + lambda_s);
                std_c = 1/sqrt(lambda + lambda_s);
                term = zeros(size(x));
                
                for jj = 1:compute.setsize
                    term(jj,:) = ((sum(x)-x(jj,:))*lambda).^2/((compute.setsize-1)*lambda + lambda_s)/2 - (sum(x.^2) - x(jj,:).^2)*lambda/2;
                end
                obs_response = sum((1-2*normcdf(0,x_c, std_c)).*exp(-x.^2/(sigma(kk)^2+std_s^2)/2).*exp(term));
                obs_response(obs_response>0)=1; obs_response(obs_response<=0)=-1;
                obs_response_opt(kk,:) = obs_response;
            end
            for mm = model_idx
                tic
                pred2opt = zeros(size(sigma));
                for kk = 1:length(sigma)
                    lambda = 1/sigma(kk)^2;
                    x = squeeze(xMat(:,:,kk));
                     % decision rule
                    obs_response = compute.decision_rule(mm,x,lambda,length(xMat));
                    obs_response(obs_response>0) = 1; obs_response(obs_response<=0) = -1;
                    
                    pred2opt(kk) = mean(obs_response==obs_response_opt(kk,:));
                end
                dirname = compute.model_names{mm};
                if ~exist([compute.dirs dirname '/pred2opt/'], 'dir');
                        mkdir([compute.dirs dirname '/pred2opt/']);
                end
                save([compute.dirs dirname '/pred2opt/pred2opt.mat'  ], 'pred2opt','sigma');
                t = toc;
                ETL = (model_idx(end) - mm)*t;
                disp(['ETL=' num2str(ETL) 's']);
            end
        end
        function prediction_calculation(model_idx)
            sigma=0;
            load('measurements.mat');
            for mm = model_idx
                tic
                predictionMat = zeros(length(sigma),length(xMat));
                for kk = 1:length(sigma)
                    lambda = 1/sigma(kk)^2;
                    x = squeeze(xMat(:,:,kk));
                     % decision rule
                    obs_response = compute.decision_rule(mm,x,lambda,length(xMat));
                    obs_response(obs_response>0) = 1; obs_response(obs_response<=0) = -1;
                    predictionMat(kk,:) = obs_response;
                end
                dirname = compute.model_names{mm};
                if ~exist('similarity','dir');
                    mkdir('similarity');
                end
                save(['similarity/' dirname],'predictionMat');
                toc
            end
        end
        function prediction_dissimilarity
            similarity = ones(25,25);
            for ii = 1:24
                dirname = compute.model_names{ii};
                load(['similarity/' dirname '.mat']);
                prediction_template = predictionMat;
                for jj = ii:25
                    dirname2 = compute.model_names{jj};
                    load(['similarity/' dirname2 '.mat']);
                    temp = predictionMat==prediction_template;
                    similarity(ii,jj) = mean(temp(:));
                    similarity(jj,ii) = similarity(ii,jj);
                end
                
            end
            distMat = 1-similarity;
            save('similarity/dissimilarity.mat','distMat','similarity')
        end
        function fake_data_generation(model_idx,Runs,flex_s)
            % generate fake data, does not entirely work with flex_s
            % 15-07-08 SS
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            if flex_s
                common_dir = compute.dirs_fakedata2;
            else
                common_dir = compute.dirs_fakedata1;
            end
            for mm = model_idx
                if flex_s
                    if sum(compute.invalid_idx==mm)
                        disp('Cannot generate fake data with this model');
                        continue
                    end
                end
                for nn = Runs
                    lambda = min(compute.lambda_fake) + range(compute.lambda_fake)*rand;
                    guess_rate = min(compute.guess_fake) + range(compute.guess_fake)*rand;
                    if flex_s
                        s_std = min(compute.std_s_vec) + range(compute.std_s_vec)*rand;
                    end
                    filename = [compute.model_names{mm} '_' num2str(nn)];
                    
                    load('fake_data_stimuli');
                    trial_num = length(stimuli_fake);
                    stimuli_fake = [stimuli_fake(1,:);repmat(stimuli_fake(2,:),compute.setsize-1,1)];
                    sigma = 1/sqrt(lambda);
                    guess_num = binornd(trial_num, guess_rate);

                    %% generate fake data
                    % internal representation of the stimuli are the stimuli plus the noise
                    x = stimuli_fake + normrnd(zeros(size(stimuli_fake)), sigma);
                    
                    % decision rule
                    if flex_s
                        obs_response = compute.decision_rule(mm,x,lambda,length(stimuli_fake),s_std);
                    else
                        obs_response = compute.decision_rule(mm,x,lambda,length(stimuli_fake));
                    end
                    
                    if mm==25
                        testor = rand(size(obs_response));
                        p_right = obs_response;
                        obs_response(testor<p_right)=1;
                        obs_response(testor>p_right)=-1;
                        obs_response(testor==p_right)=0;
                    else
                        obs_response = sign(obs_response);
                    end
                    temp = obs_response(obs_response==0);
                    if ~isempty(temp)
                        obs_response(obs_response==0)= 2*randi(2,size(temp))-3;
                    end
                    idx = randi(trial_num-guess_num);
                    obs_response(idx:idx+guess_num-1) = 2*randi(2,size(obs_response(idx:idx+guess_num-1)))-3;

                    % store fake data
                    data = [stimuli_fake' obs_response'];
                    if ~exist(common_dir,'dir')
                        mkdir(common_dir);
                    end
                    if flex_s
                        save([common_dir filename], 'data','lambda','guess_rate','s_std');
                    else
                        save([common_dir filename], 'data','lambda','guess_rate');
                    end
                end
            end
        end
        function fake_data_generation_real_pars(model_idx,flex_s)
            % generate fake data, with real parameters from the subjects, does not entirely work with flex_s
            % 15-07-09 SS
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            if flex_s
                common_dir = compute.dirs_fakedata2;
            else
                common_dir = compute.dirs_fakedata1;
            end
            Runs = 1:9; % to match the real subjects
            for mm = model_idx
                if flex_s
                    if sum(compute.invalid_idx==mm)
                        disp('Cannot generate fake data with this model');
                        continue
                    end
                end
                for nn = Runs
                    % load weighted parameters from real subjects
                    [lambda, guess_rate] = utils.generate_fake_pars(nn);
               
                    if flex_s
                        s_std = min(compute.std_s_vec) + range(compute.std_s_vec)*rand;
                    end
                    filename = [compute.model_names{mm} '_' num2str(nn)];
                    
                    load('fake_data_stimuli');
                    trial_num = length(stimuli_fake);
                    stimuli_fake = [stimuli_fake(1,:);repmat(stimuli_fake(2,:),compute.setsize-1,1)];
                    sigma = 1/sqrt(lambda);
                    guess_num = binornd(trial_num, guess_rate);

                    %% generate fake data
                    % internal representation of the stimuli are the stimuli plus the noise
                    x = stimuli_fake + normrnd(zeros(size(stimuli_fake)), sigma);
                    
                    % decision rule
                    if flex_s
                        obs_response = compute.decision_rule(mm,x,lambda,length(stimuli_fake),s_std);
                    else
                        obs_response = compute.decision_rule(mm,x,lambda,length(stimuli_fake));
                    end
                    
                    if mm==28 % sampling model
                        testor = rand(size(obs_response));
                        p_right = obs_response;
                        obs_response(testor<p_right)=1;
                        obs_response(testor>p_right)=-1;
                        obs_response(testor==p_right)=0;
                    else
                        obs_response = sign(obs_response);
                    end
                    temp = obs_response(obs_response==0);
                    if ~isempty(temp)
                        obs_response(obs_response==0)= 2*randi(2,size(temp))-3;
                    end
                    idx = randi(trial_num-guess_num);
                    obs_response(idx:idx+guess_num-1) = 2*randi(2,size(obs_response(idx:idx+guess_num-1)))-3;

                    % store fake data
                    data = [stimuli_fake' obs_response'];
                    if ~exist(common_dir,'dir')
                        mkdir(common_dir);
                    end
                    if flex_s
                        save([common_dir filename], 'data','lambda','guess_rate','s_std');
                    else
                        save([common_dir filename], 'data','lambda','guess_rate');
                    end
                end
            end
        end
        function likelihood_calculation(subj_idx,run_idx,model_idx,flex_s)
            % calculate likelihood of giving a certain answer for fake
            % data, based on the pre computed big prediction table
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            if flex_s
                common_dir = compute.dirs_fake2;
                common_dir2 = compute.dirs_fakedata2;
            else
                common_dir = compute.dirs_fake1;
                common_dir2 = compute.dirs_fakedata1;
            end
            trial_num_sim = 1000;
            for mm = subj_idx
                for jj = run_idx
                    for kk = model_idx
                        dirname = compute.model_names{kk};
                        subjid = [compute.model_names{mm} '_' num2str(jj)];
                        load([common_dir dirname '/tables.mat']);
                        load([common_dir2 subjid]);
                        response = data(:,5);
                        lambda = compute.lambda_vec;
                        guess = compute.guess_vec;
                        if flex_s
                            s_std = compute.std_s_vec;
                            log_CPG_prediction = zeros(length(lambda),length(s_std),length(guess));
                            CP_prediction_r(response==-1,:) = 1 - CP_prediction_r(response==-1,:);
                        else
                            log_CPG_prediction = zeros(length(lambda), length(guess));
                            CP_prediction_r(response==-1,:,:) = 1 - CP_prediction_r(response==-1,:,:);
                        end
                        for ii = 1:length(guess)
                            CPG_prediction_r = CP_prediction_r*(1-guess(ii)) + 0.5*guess(ii);
                            CPG_prediction_r(CPG_prediction_r==0) = 1/trial_num_sim/10;
                            CPG_prediction_r(CPG_prediction_r==1) = 1 - 1/trial_num_sim/10;
                            if flex_s
                                log_CPG_prediction(:,:,ii) = squeeze(sum(log(CPG_prediction_r),1));
                            else
                                log_CPG_prediction(:,ii) = squeeze(sum(log(CPG_prediction_r),1));
                            end
                        end
                        
                        if ~exist([common_dir dirname '/tables'], 'dir');
                            mkdir([common_dir dirname '/tables']);
                        end
                        save([common_dir dirname '/tables/' subjid], 'log_CPG_prediction');
                     end
                end
            end
        end
        function table_bin(model_idx)
            % compute a table with prediction of reporting right for each
            % bin, the size of the returned table is
            % length(bins),length(bins),length(lambda),length(guess)
            % 2015-08-20 SS
            for kk = model_idx
                dirname = compute.model_names{kk};
                load([compute.dirs dirname '/tables_grid.mat']);
                lambda = compute.lambda_vec;
                guess = compute.guess_vec;

                guessMat = repmat(guess',[1,length(compute.sT_vec), length(compute.sD_vec), length(lambda)]);
                guessMat = permute(guessMat,[2,3,4,1]);
                CP_prediction_r_expand = repmat(CP_prediction_r,[1,1,1,length(guess)]);
                CPG_prediction_r = CP_prediction_r_expand.*(1-guessMat)+0.5*guessMat;
                
                bins = utils.generate_bins(9);
                idx1 = interp1(bins,1:length(bins),compute.sT_vec, 'nearest','extrap');
                idx2 = interp1(bins,1:length(bins),compute.sD_vec, 'nearest','extrap');
                CPG_prediction_mean = zeros(length(bins),length(bins),length(lambda),length(guess));
                % get the mean predictions within bins
                for ii = 1:length(bins)
                    for jj = 1:length(bins)
                        prediction_rel = CPG_prediction_r(idx1==ii,idx2==jj,:,:);
                        pStim = normpdf(compute.sT_vec(idx1==ii)',0,compute.s_std)*normpdf(compute.sD_vec(idx2==jj),0,compute.s_std);
                        pStimSum = sum(pStim(:));
                        pStimMat = repmat(pStim, [1,1,length(lambda),length(guess)]);
                        prediction_mean = squeeze(sum(squeeze(sum(pStimMat.*prediction_rel))))/pStimSum;
                        CPG_prediction_mean(ii,jj,:,:) = prediction_mean;
                    end
                end
                CPG_prediction_mean(CPG_prediction_mean>=1) = 1-1/10000;
                CPG_prediction_mean(CPG_prediction_mean==0) = 1/10000;
                
                save([compute.dirs dirname '/tables_bins'], 'CPG_prediction_mean')
            end
                
        end
        function evi_bin_grid(subj_idx, model_idx)
            
            for mm = subj_idx
                subjid = compute.subjids{mm};
                load([compute.dirs '/data_bin/' subjid '.mat'])
                for kk = model_idx
                    dirname = compute.model_names{kk};
                    load([compute.dirs dirname '/tables_bins'])
                    cntMat = repmat(cnt, [1,1,length(compute.lambda_vec),length(compute.guess_vec)]);
                    cnt_rMat = repmat(cnt_r, [1,1,length(compute.lambda_vec),length(compute.guess_vec)]);
                    cnt_lMat = cntMat - cnt_rMat;
                    bmcMat = log(CPG_prediction_mean).*cnt_rMat + log(1-CPG_prediction_mean).*cnt_lMat;
                    % sum over different parameters
                    
                    bmcMat = squeeze(sum(squeeze(sum(bmcMat))));
                    L_max = max(bmcMat(:));
                    L = bmcMat - L_max;
                    bmc = L_max + log(sum(exp(L(:)))) - log(numel(bmcMat));
                    bic = L_max - log(2000);
                    aic = L_max - 2;
                    
                    if ~exist([compute.dirs dirname '/evi_bin_grid'], 'dir');
                        mkdir([compute.dirs dirname '/evi_bin_grid']);
                    end
                    save([compute.dirs dirname '/evi_bin_grid/' subjid], 'bmc','bic','aic','bmcMat');
                end
            end

        end
        function mle_bin_grid(subj_idx, model_idx, cross)
            % compute maximum log likelihood and its 95% confidence interval of a model
            % cross specifies whether using cross validation
            % 0: no cross validation, use the whole data set
            % 1: with cross validation, use the first half of data to get
            % the best parameters and calculate mle with the second half of
            % data
            
            % SS 2015-08-24
            
            if ~exist('cross','var')
                cross = 0;
            end
            
            if cross == 0
                data_dir = 'data_bin_0';
                half = 0;
            else
                data_dir = 'data_bin_1';
                half = 2;
            end
            
            for mm = subj_idx
                
                subjid = compute.subjids{mm};
                load([compute.dirs data_dir '/' subjid '.mat'])
                
                % check whether entropy of this subject has been computed.
                if ~exist(['real_data_results/entropy_est_' num2str(half) '/' subjid '_9.mat'],'file')
                    error('Compute the entropy first!')
                end
                
                load(['real_data_results/entropy_est_' num2str(half) '/' subjid '_9.mat']);
                
                for kk = model_idx
                    dirname = compute.model_names{kk};
                    load([compute.dirs dirname '/tables_bins'])
                    cntMat = repmat(cnt, [1,1,length(compute.lambda_vec),length(compute.guess_vec)]);
                    cnt_rMat = repmat(cnt_r, [1,1,length(compute.lambda_vec),length(compute.guess_vec)]);
                    cnt_lMat = cntMat - cnt_rMat;
                    predMat = log(CPG_prediction_mean).*cnt_rMat + log(1-CPG_prediction_mean).*cnt_lMat;
                    
                    % pick out the best parameters
                    
                    predMat2 = squeeze(sum(squeeze(sum(predMat))));
                    [~,I] = max(predMat2(:));
                    [idx1,idx2] = ind2sub([length(compute.lambda_vec),length(compute.guess_vec)],I);
                    lambda_hat = compute.lambda_vec(idx1);
                    guess_hat = compute.lambda_vec(idx2);
                    CPG_prediction_hat = CPG_prediction_mean(:,:,idx1,idx2);
                        
                    if cross == 1 % with cross validation, recompute mle with second half of data and best pars got from the first half
                        load([compute.dirs '/data_bin_2/' subjid '.mat'])
                        cnt_l = cnt - cnt_r;
                        mleMat = log(CPG_prediction_hat).*cnt_r + log(1-CPG_prediction_hat).*cnt_l;
                    end
                    
                    mle_mean = sum(mleMat(:));    
                    varMat = (cnt_r+1).*(cnt-cnt_r+1)./(cnt+2).^2./(cnt+3);
                    mle_varMat = cnt.^2.*(log(CPG_prediction_hat)-log(1-CPG_prediction_hat)).^2.*varMat ;
                    mle_var = sum(mle_varMat(:));
                    mle_err = norminv(0.975)*sqrt(mle_var);
                    mle_range = [mle_mean - mle_err, mle_mean + mle_err];
                    dkl = entropyvalue - mle_mean;
                    % compute the p value of how distinguishable llmax with
                    % the entropy
                    p = 1-normcdf((entropyvalue-mle_mean)/sqrt(mle_var));
                    
                    if cross == 0
                        save_dir = '/mle_bin_grid';
                    else
                        save_dir = '/mle_bin_grid_cross';
                    end
                    
                    if ~exist([compute.dirs dirname save_dir], 'dir');
                        mkdir([compute.dirs dirname save_dir]);
                    end
                    save([compute.dirs dirname save_dir '/' subjid], 'mleMat','mle_mean','mle_var','mle_range','mle_err','p','dkl','CPG_prediction_hat','lambda_hat','guess_hat');
                end
            end
            
        end
        function simulation_homo(sigma)
            % a simulation to compare the model agreement of max rule and
            % opt rule in homogeneous experiment
            % SS 03-04-2014
            load('stimuli.mat')
            stimuli = [sT; repmat(sD,compute.setsize-1,1)];
            xMat = zeros(size(stimuli,1), size(stimuli,2),length(sigma));
            for ii = 1:length(sigma);
                xMat(:,:,ii) = normrnd(stimuli, sigma(ii));
            end
            
            pred2opt = zeros(size(sigma));
            for ii = 1:length(sigma)
                x = squeeze(xMat(:,:,ii));
                % calculate prediction of the optimal model
                lambda = 1/sigma(ii)^2;
                maxx = max(x.^2);
                f = exp(bsxfun(@minus,x.^2,maxx)*lambda/2);
                x_c = lambda*x/(lambda + 1/compute.s_std^2);
                lambda_c = lambda + 1/compute.s_std^2;
                term1 = squeeze(sum(f.*(1-normcdf(0,x_c,1/sqrt(lambda_c)))));
                term2 = squeeze(sum(f.*normcdf(0,x_c,1/sqrt(lambda_c))));
                obs_response_opt = sign(term1 - term2);
                % calculate prediction of max model
                [~,idx] = max(abs(x));
                idx = sub2ind(size(x), idx, 1:length(x));
                obs_response_max = x(idx);
                obs_response_max = sign(obs_response_max);
                pred2opt(ii) = mean(obs_response_max==obs_response_opt);
            end
            save('homo_pred2opt', 'pred2opt','sigma');
        end
        function [sigmaEst, sigma_sEst]=simulation_learning(nTrials,sigma)
            % simulate a process of learning, check whether estimated
            % sigma, and sigma_s will converge to the true value across
            % trials.
            % SS 03-05-2014
            sigmaVec = linspace(1,20,100);
            sigma_sVec = linspace(1,20,100);
            sigmaMat_temp= repmat(sigmaVec',1,length(sigma_sVec));
            sigma_sMat_temp = repmat(sigma_sVec,length(sigmaVec),1);
            sigmaMat = repmat(permute(sigmaMat_temp,[3,1,2]),[compute.setsize,1,1]);
            sigma_sMat = repmat(permute(sigma_sMat_temp,[3,1,2]),[compute.setsize,1,1]);
            lambdaMat = 1/sigmaMat.^2;
            lambda_sMat = 1/sigma_sMat.^2;
            LLCum = zeros(size(sigmaMat_temp));
            sigmaEst = zeros(1,nTrials);
            sigma_sEst = zeros(1,nTrials);
            for ii = 1:nTrials
                sT = normrnd(0,compute.s_std);
                sD = normrnd(0,compute.s_std);
                stimulus = [sT,sD,sD,sD]';
                x = normrnd(stimulus,sigma);
                xMat = repmat(x,[1,length(sigmaVec),length(sigma_sVec)]);
                term_x = zeros(size(xMat));
                idx_temp = 1:compute.setsize;
                for jj = 1:compute.setsize
                    temp1 = squeeze(xMat(idx_temp==jj,:,:));
                    temp2 = squeeze(xMat(idx_temp~=jj,:,:));
                    term_x(jj,:,:) = temp1.^2./(sigmaMat_temp.^2 + sigma_sMat_temp.^2) + squeeze(mean(temp2)).^2./(sigmaMat_temp.^2/(compute.setsize-1) + sigma_sMat_temp.^2) ...
                        + squeeze(var(temp2,1)).*(compute.setsize-1)./sigmaMat_temp.^2;
                end
                
                term = lambdaMat.^2.*lambda_sMat.*sqrt(1./(lambdaMat+lambda_sMat)).*sqrt(1./((compute.setsize-1)*lambdaMat+lambda_sMat));
                if sT>0
                    p = squeeze(sum((1+erf(xMat.*lambdaMat./sqrt(2*(lambdaMat+lambda_sMat)))).*exp(-term_x/2).*term));
                else
                    p = squeeze(sum((1-erf(xMat.*lambdaMat./sqrt(2*(lambdaMat+lambda_sMat)))).*exp(-term_x/2).*term));
                end
                p(p==0) = 1e-100;
                p(p==1) = 1-1e-100;
                LLCum = LLCum + log(p);
                [~,idx] = max(LLCum(:));
                [idx1,idx2] = ind2sub(size(LLCum),idx);
                sigmaEst(ii) = sigmaVec(idx1);
                sigma_sEst(ii) = sigma_sVec(idx2);
%                 figure(101);
%                 imagesc(LLCum);
            end
%             figure; plot(1:nTrials,sigmaEst,1:nTrials,sigma_sEst);
            
        end
        function [sigmaEst, sigma_sEst]=learning(stimuli,sigma)
            % given a series of simulus, return the learned sigma and
            % sigma_s after each trial
            % SS 03-28-2014
            nTrials = length(stimuli);
            sigmaVec = linspace(1,20,100);
            sigma_sVec = linspace(1,20,100);
            sigmaMat_temp= repmat(sigmaVec',1,length(sigma_sVec));
            sigma_sMat_temp = repmat(sigma_sVec,length(sigmaVec),1);
            sigmaMat = repmat(permute(sigmaMat_temp,[3,1,2]),[compute.setsize,1,1]);
            sigma_sMat = repmat(permute(sigma_sMat_temp,[3,1,2]),[compute.setsize,1,1]);
            lambdaMat = 1/sigmaMat.^2;
            lambda_sMat = 1/sigma_sMat.^2;
            LLCum = zeros(size(sigmaMat_temp));
            sigmaEst = zeros(1,nTrials);
            sigma_sEst = zeros(1,nTrials);
            for ii = 1:nTrials
                sT = stimuli(ii,1);
                sD = stimuli(ii,2);
                stimulus = [sT,sD,sD,sD]';
                x = normrnd(stimulus,sigma);
                xMat = repmat(x,[1,length(sigmaVec),length(sigma_sVec)]);
                term_x = zeros(size(xMat));
                idx_temp = 1:compute.setsize;
                for jj = 1:compute.setsize
                    temp1 = squeeze(xMat(idx_temp==jj,:,:));
                    temp2 = squeeze(xMat(idx_temp~=jj,:,:));
                    term_x(jj,:,:) = temp1.^2./(sigmaMat_temp.^2 + sigma_sMat_temp.^2) + squeeze(mean(temp2)).^2./(sigmaMat_temp.^2/(compute.setsize-1) + sigma_sMat_temp.^2) ...
                        + squeeze(var(temp2,1)).*(compute.setsize-1)./sigmaMat_temp.^2;
                end
                
                term = lambdaMat.^2.*lambda_sMat.*sqrt(1./(lambdaMat+lambda_sMat)).*sqrt(1./((compute.setsize-1)*lambdaMat+lambda_sMat));
                if sT>0
                    p = squeeze(sum((1+erf(xMat.*lambdaMat./sqrt(2*(lambdaMat+lambda_sMat)))).*exp(-term_x/2).*term));
                else
                    p = squeeze(sum((1-erf(xMat.*lambdaMat./sqrt(2*(lambdaMat+lambda_sMat)))).*exp(-term_x/2).*term));
                end
                p(p==0) = 1e-100;
                p(p==1) = 1-1e-100;
                LLCum = LLCum + log(p);
                [~,idx] = max(LLCum(:));
                [idx1,idx2] = ind2sub(size(LLCum),idx);
                sigmaEst(ii) = sigmaVec(idx1);
                sigma_sEst(ii) = sigma_sVec(idx2);
%                 figure(101);
%                 imagesc(LLCum);
            end
%             figure; plot(1:nTrials,sigmaEst,1:nTrials,sigma_sEst);
            
        end
        function simulation_learning_ntimes(nTrials,sigma,n)
            % run simulation for multiple times, show the mean and s.e.m of
            % the covergence curve.
            % SS 03-06-2014
            sigEstMat = zeros(n,nTrials);
            sigsEstMat = zeros(n,nTrials);
            for ii = 1:n
                [sigEstMat(ii,:),sigsEstMat(ii,:)] = compute.simulation_learning(nTrials,sigma);
            end
            save('learning', 'sigEstMat', 'sigsEstMat', 'sigma', 'nTrials');
        end
        function simulation_learning_fix_sigs(nTrials,sigma)
            % with a fixed sigma_s, test whether sigma could be learned across trials
            % SS 03-06-2014
            sigmaVec = linspace(1,20,100);
            sigma_s = compute.s_std;
            sigmaMat = repmat(sigmaVec,compute.setsize,1);
            lambdaMat = 1./sigmaMat.^2;
            lambda_s = 1/sigma_s^2;
            LLCum = zeros(size(sigmaVec));
            sigmaEst = zeros(1,nTrials);
            for ii = 1:nTrials
                sT = normrnd(0,compute.s_std);
                sD = normrnd(0,compute.s_std);
                stimulus = [sT,sD,sD,sD]';
                x = normrnd(stimulus,sigma);
                xMat = repmat(x,1,length(sigmaVec));
                term_x = zeros(size(xMat));
                idx_temp = 1:compute.setsize;
                for jj = 1:compute.setsize
                    temp1 = squeeze(xMat(idx_temp==jj,:));
                    temp2 = squeeze(xMat(idx_temp~=jj,:));
                    term_x(jj,:) = temp1.^2./(sigmaVec.^2 + sigma_s^2) + squeeze(mean(temp2)).^2./(sigmaVec.^2/(compute.setsize-1) + sigma_s^2) ...
                        + squeeze(var(temp2)).*(compute.setsize-1)./sigmaVec.^2;
                end
                term = lambdaMat.^2.*lambda_s.*sqrt(1./(lambdaMat+lambda_s)).*sqrt(1./((compute.setsize-1)*lambdaMat+lambda_s));
                x_c = xMat.*lambdaMat./(lambdaMat + lambda_s);
                std_c = 1./sqrt(lambdaMat + lambda_s);
                
                if sT>=0
                    p = squeeze(sum((1-normcdf(0,x_c,std_c)).*exp(-term_x/2).*term));
                else
                    p = squeeze(sum(normcdf(0,x_c,std_c).*exp(-term_x/2).*term));
                end
                p(p==0) = 1e-100;
                p(p==1) = 1-1e-100;
                LLCum = LLCum + log(p);
                [~,idx] = max(LLCum(:));
                sigmaEst(ii) = sigmaVec(idx);
%               figure(101);
%                 imagesc(LLCum);
            end
            figure; plot(1:nTrials,sigmaEst);
        end
        function simulation_learning_fix_sig(nTrials,sigma)
            % with a fixed sigma, test whether sigma_s could be learned across trials
            % SS 03-17-2014
            sigma_sVec = linspace(1,20,100);
            sigma_sMat = repmat(sigma_sVec,compute.setsize,1);
            lambda = 1/sigma^2;
            lambda_sMat = 1./sigma_sMat.^2;
            LLCum = zeros(size(sigma_sVec));
            sigma_sEst = zeros(1,nTrials);
            for ii = 1:nTrials
                sT = normrnd(0,compute.s_std);
                sD = normrnd(0,compute.s_std);
                stimulus = [sT,sD,sD,sD]';
                x = normrnd(stimulus,sigma);
                xMat = repmat(x,1,length(sigma_sVec));
                term_x = zeros(size(xMat));
                idx_temp = 1:compute.setsize;
                for jj = 1:compute.setsize
                    temp1 = squeeze(xMat(idx_temp==jj,:));
                    temp2 = squeeze(xMat(idx_temp~=jj,:));
                    term_x(jj,:) = temp1.^2./(sigma^2 + sigma_sVec.^2) + squeeze(mean(temp2)).^2./(sigma^2/(compute.setsize-1) + sigma_sVec.^2) ...
                        + squeeze(var(temp2)).*(compute.setsize-1)./sigma^2;
                end
                term = lambda^2.*lambda_sMat.*sqrt(1./(lambda+lambda_sMat)).*sqrt(1./((compute.setsize-1)*lambda+lambda_sMat));
                x_c = xMat.*lambda./(lambda + lambda_sMat);
                std_c = 1./sqrt(lambda + lambda_sMat);
 

                if sT>=0
                    p = squeeze(sum((1-normcdf(0,x_c,std_c)).*exp(-term_x/2).*term));
                else
                    p = squeeze(sum(normcdf(0,x_c,std_c).*exp(-term_x/2).*term));
                end
                p(p==0) = 1e-100;
                p(p==1) = 1-1e-100;
                LLCum = LLCum + log(p);
                [~,idx] = max(LLCum(:));
                sigma_sEst(ii) = sigma_sVec(idx);
%               figure(101);
%                 imagesc(LLCum);
            end
            figure; plot(1:nTrials,sigma_sEst);
            sigma_sEst(end)
        end
        function simulation_learning_one_item(nTrials,sigma)
            % with a fixed sigma_s, test whether sigma could be learned
            % across trials, for simple one item discrimination experiment
            % SS 03-17-2014
            sigmaVec = linspace(1,20,100);
            sigma_s = compute.s_std;
            lambdaMat = 1./sigmaVec.^2;
            lambda_s = 1/sigma_s^2;
            LLCum = zeros(size(sigmaVec));
            sigmaEst = zeros(1,nTrials);
            for ii = 1:nTrials
                stimulus = normrnd(0,compute.s_std);
                x = normrnd(stimulus,sigma);
                xMat = repmat(x,1,length(sigmaVec));                        
                term_x = xMat.^2./(sigmaVec.^2 + sigma_s^2);
                term = sqrt(lambdaMat.*lambda_s).*sqrt(1./(lambdaMat+lambda_s));
                x_c = xMat.*lambdaMat./(lambdaMat + lambda_s);
                std_c = 1./sqrt(lambdaMat + lambda_s);
                
                if stimulus>=0
                    p = (1-normcdf(0,x_c,std_c)).*exp(-term_x/2).*term;
                else
                    p = normcdf(0,x_c,std_c).*exp(-term_x/2).*term;
                end
                p(p==0) = 1e-100;
                p(p==1) = 1-1e-100;
                LLCum = LLCum + log(p);
                [~,idx] = max(LLCum(:));
                sigmaEst(ii) = sigmaVec(idx);
%               figure(101);
%                 imagesc(LLCum);
            end
            figure; plot(1:nTrials,sigmaEst); 
        end
        function simulation_learning_homo(nTrials,sigma)
            % with a fixed sigma_s, test whether sigma could be learned
            % across trials, for homogeneous experiment (3 vertical)
            % SS 03-18-2014
            sigmaVec = linspace(1,20,100);
            sigma_s = compute.s_std;
            sigmaMat = repmat(sigmaVec,compute.setsize,1);
            lambdaMat = 1./sigmaMat.^2;
            lambda_s = 1/sigma_s^2;
            LLCum = zeros(size(sigmaVec));
            sigmaEst = zeros(1,nTrials);
            for ii = 1:nTrials
                sT = normrnd(0,compute.s_std);
                sD = 0;
                stimulus = [sT,sD,sD,sD]';
                x = normrnd(stimulus,sigma);
                xMat = repmat(x,1,length(sigmaVec));
                term_x = zeros(size(xMat));
                idx_temp = 1:compute.setsize;
                for jj = 1:compute.setsize
                    temp1 = squeeze(xMat(idx_temp==jj,:));
                    temp2 = squeeze(xMat(idx_temp~=jj,:));
                    term_x(jj,:) = temp1.^2./(sigmaVec.^2 + sigma_s^2) + squeeze(sum(temp2.^2))./sigmaVec.^2;
                end
                term = lambdaMat.^2.*lambda_s.*sqrt(1./(lambdaMat+lambda_s));
                x_c = xMat.*lambdaMat./(lambdaMat + lambda_s);
                std_c = 1./sqrt(lambdaMat + lambda_s);
                
                if sT>=0
                    p = squeeze(sum((1-normcdf(0,x_c,std_c)).*exp(-term_x/2).*term));
                else
                    p = squeeze(sum(normcdf(0,x_c,std_c).*exp(-term_x/2).*term));
                end
                p(p==0) = 1e-100;
                p(p==1) = 1-1e-100;
                LLCum = LLCum + log(p);
                [~,idx] = max(LLCum(:));
                sigmaEst(ii) = sigmaVec(idx);
%               figure(101);
%                 imagesc(LLCum);
            end
            figure; plot(1:nTrials,sigmaEst);
        end
        function predictions_learning(subj_idx,session, flex_s)
            if ~exist('flex_s','var')
                flex_s = 0;
            end
            if flex_s
                common_dir = compute.dirs3;
            else
                common_dir = compute.dirs;
            end
            nn=1;           
            for mm = subj_idx
                % calculate fit predictions for models and p_right for real data
                dirname = compute.model_names{1};
                subjid = compute.subjids{mm};
                %% load data
                load([common_dir dirname '/fit_pars/' subjid]);
                [stimuli,~,performance] = utils.readdata2(subjid,session);
                trial_num = length(stimuli);

                trial_num_sim = 1000;
                [sigmaEst, sigma_sEst] = compute.learning(stimuli,1/sqrt(CPG_lambda_hat));
                %% pars settings
                stimuli = [stimuli(:,1) repmat(stimuli(:,2),1,compute.setsize-1)];
                sT = stimuli(:,1)';
                %% CP & CPG models
                % loop over real trials, calculate the prediction for each stimulus
                for i = 1:trial_num
                    noiseMat_CPG = normrnd(0,repmat(1/sqrt(CPG_lambda_hat),compute.setsize,trial_num_sim));

                    %% CPG
                    % compute predictions
                    x = repmat(stimuli(i,:),trial_num_sim,1);
                    x = x' + noiseMat_CPG;
                    lambda_temp = 1/sigmaEst(i)^2;
                    std_s = sigma_sEst(i);
                    % decision rule
                    if flex_s
                        obs_response = compute.decision_rule(nn,x,lambda_temp,trial_num_sim,std_s);
                    else
                        obs_response = compute.decision_rule(nn,x,lambda_temp,trial_num_sim);
                    end
                        % the prediction of the probability of reporting right
               
                    CPG_prediction(i) = (sum(obs_response>0) + .5*sum(obs_response==0))/trial_num_sim;


                end
                CPG_prediction = CPG_prediction*(1-CPG_guess_hat) + .5*CPG_guess_hat;
                CPG_prediction(sT<0) = 1 - CPG_prediction(sT<0);

                %% save result
                if ~exist([common_dir dirname '/fit_results_learning/'], 'dir')
                    mkdir([common_dir dirname '/fit_results_learning/']);
                end
                save([common_dir dirname '/fit_results_learning/' subjid '_' num2str(session)], 'CPG_prediction','performance');
            end
        end
        function calculate_entropy(subj_idx,nBins)
            for mm = subj_idx
                subjid = compute.subjids{mm};
                
                [stimuli,response] = utils.readdata(subjid);
                target_stimuli = stimuli(:,1);
                dist_stimuli = stimuli(:,2);
                     
                target_stimuli_right = target_stimuli(response==1);
                dist_stimuli_right = dist_stimuli(response==1);
                
                for nn = nBins
                    bins = utils.generate_bins(nn);
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
                    %% calculate the real probability of reporting right
                    p_right = cnt_r./cnt;
                    cntVec = cnt;

                    entropyMat = p_right.*log(p_right)+(1-p_right).*log(1-p_right);
                    entropyMat(p_right==0) = 0;
                    entropyMat(p_right==1) = 0;
                    entropyvalue = sum(entropyMat(:).*cntVec(:));
                    save_dir = 'real_data_results/entropy/';
                    if ~exist(save_dir,'dir')
                        mkdir(save_dir)
                    end
                    save([save_dir subjid '_' num2str(nn)],'entropyMat','entropyvalue');
                end
                
             end
        end
        function calculate_entropy_est(subj_idx,nBins,half)
            % Apply estimation of entropy from Grassberger, 2008
            % half specifies which part of data to use
            % 0: whole data set
            % 1: the half with odd indices
            % 2: the half with even indices
            % 2015-08-20 SS
            % modified 2015-09-14 SS
            for mm = subj_idx
                subjid = compute.subjids{mm};
                
                [stimuli,response] = utils.readdata(subjid,half);
                                  
                for nn = nBins
                    bins = utils.generate_bins(nn);
                    [cnt,cnt_r] = utils.compute_cnts(stimuli,response,bins);

                    [entropyvalue,entropyMat] = utils.compute_G_entropy(cnt,cnt_r);
                    save_dir = ['real_data_results/entropy_est_' num2str(half) '/'];
                    if ~exist(save_dir,'dir')
                        mkdir(save_dir)
                    end
                    save([save_dir subjid '_' num2str(nn)],'entropyMat','entropyvalue');
                end
                
             end
        end
        function calculate_Dkl_err_bootstrap(subj_idx,model_idx)
            % bootstrapping method to get the estimation error for KL
            % divergence, now only works for the cross-validated data
            nRep = 10000;
            for ii = 1:length(subj_idx)
                subjid = compute.subjids{subj_idx(ii)};
                load([compute.dirs 'data_bin_2/' subjid '.mat'],'p_right','cnt_r','cnt')
                for jj = 1:length(model_idx)
                    dirname = compute.model_names{model_idx(jj)};
                    load([compute.dirs dirname '/mle_bin_grid_cross/' subjid],'CPG_prediction_hat');
                    dkl = zeros(1,nRep);
                    for iRep = 1:nRep
                        
                        cnt_r_temp = binornd(cnt, p_right);
                        mle = utils.compute_cross_entropy(cnt,cnt_r_temp,CPG_prediction_hat);
                        entropy = utils.compute_G_entropy(cnt,cnt_r_temp);
                        dkl(iRep) = entropy - mle;
                    end
                    figure; hist(dkl);
                    pvalue = sum(dkl<0)/nRep
                    dkl_ci_95 = [prctile(dkl,2.5),prctile(dkl,97.5)]
                    save(['real_data_results/entropy_est_2/' subjid '_9_bootstrap'],'pvalue','dkl_ci_95');
                end
            end
        end
        
    end
end
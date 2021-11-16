function dec_rule = calculate_optimaldecisioncriteria(pars)

persistent optimal_criterion
persistent saved_pars
persistent k_range
persistent J_lin
persistent highest_J

Jbar_high = pars(1);
Jbar_low = pars(2);
tau = pars(3);

nItems = 4;
nLowVec = nItems:-1:0;
nRelConds = length(nLowVec);

if isempty(saved_pars) || any(saved_pars ~=pars)
    
    % make CDF for interpolating J to Kappa
    if isempty(k_range)
        tempp = load('cdf_table.mat');
        % K_interp = tempp.K_interp;
        % cdf = tempp.cdf;
        k_range = tempp.k_range;
        J_lin = tempp.J_lin;
        highest_J = tempp.highest_J;
    end
    
    %     figure;
    
    nSim = 2*1e5; % half of this will make up d0 and d1
    nReps = 1;
    optimal_criterion = zeros(nReps,nRelConds);
    
    for icond = 1:nRelConds
        nLow = nLowVec(icond);
        
        % mean precision, Jbar
        Jbar_mat = nan(nSim,nItems);
        Jbar_mat(:,1:nLow) = Jbar_low;
        Jbar_mat(:,(nLow+1):nItems) = Jbar_high;
        
        %         for irep = 1:nReps
        
        % sampled precision, J
        J_x_mat = gamrnd(Jbar_mat./tau,tau);
        J_y_mat = gamrnd(Jbar_mat./tau,tau);
        
        % set kappas too high to highest J (alternatively can resample, as
        % keshvari did)
        J_x_mat(J_x_mat > highest_J) = highest_J;
        J_y_mat(J_y_mat > highest_J) = highest_J;
        
        % convert J to kappa
        xi = 1/diff(J_lin(1:2))*J_x_mat+1;
        kappa_x_i = k_range(round(xi));
        xi = 1/diff(J_lin(1:2))*J_y_mat+1;
        kappa_y_i = k_range(round(xi));
        
        if size(kappa_x_i,2) ~= nItems
            kappa_x_i = kappa_x_i';
            kappa_y_i = kappa_y_i';
        end
        
        % generate measurement noise
        noise_x = circ_vmrnd(0,kappa_x_i);
        noise_y = circ_vmrnd(0,kappa_y_i);
        
        % get difference between noise
        delta_noise = noise_x-noise_y;
        
        Delta = zeros(nSim,4);
        deltaVec = 4*pi.*rand(1,nSim/2)-2*pi;
        idx = sub2ind([nSim, 4],1:(nSim/2),randi(4,1,nSim/2));
        Delta(idx) = deltaVec;
        
        % the term inside denominator bessel function for d_i
        Kc = bsxfun(@times,2.*kappa_x_i.*kappa_y_i,cos(bsxfun(@plus,Delta,delta_noise))); % note: it is okay to simply add the noise bc it goes through a cos!!
        Kc = sqrt(bsxfun(@plus,kappa_x_i.^2+kappa_y_i.^2,Kc)); % dims: mat_dims
        
        % d1 d0
        d_Mat = squeeze(max(kappa_x_i + kappa_y_i - Kc,[],2));
        d1 = d_Mat(1:nSim/2,:);
        d0 = d_Mat(nSim/2+1:end,:);
        
        % Bisection method homebrew:
        crit0 = min(d0(:));
        crit2 = max(d1(:));
        crit1 = crit0 + (crit2-crit0) / ( 3 + sqrt(5) ) * 2;
        f0 = sum(d1 > crit0) + sum(d0 < crit0); % hits, correct rejections
        f1 = sum(d1 > crit1) + sum(d0 < crit1);
        f2 = sum(d1 > crit2) + sum(d0 < crit2);
        
        while (crit2-crit0) > 0.001
            if abs(crit2-crit1) > abs(crit1-crit0)
                crit_new = crit1 + (crit2-crit1) / ( 3 + sqrt(5) ) * 2;
                f_new = sum(d1 > crit_new) + sum(d0 < crit_new);
                if f_new < f1
                    f2 = f_new;
                    crit2 = crit_new;
                else
                    f0 = f1;
                    crit0 = crit1;
                    f1 = f_new;
                    crit1 = crit_new;
                end
            else
                crit_new = crit1 + (crit0-crit1) / ( 3 + sqrt(5) ) * 2;
                f_new = sum(d1 > crit_new) + sum(d0 < crit_new);
                if f_new < f1
                    f0 = f_new;
                    crit0 = crit_new;
                else
                    f2 = f1;
                    crit2 = crit1;
                    f1 = f_new;
                    crit1 = crit_new;
                end
            end
            
        end
        %             subplot(2,3,6);
        %             histogram(Delta)
        
        %                     subplot(2,3,icond)
        %                     histogram(d0,linspace(0,10,100)); hold on; histogram(d1,linspace(0,10,100));
        %                     plot([crit1 crit1],[0 2000],'k'), pause
        %             optimal_criterion(irep,icond) = crit1;
        optimal_criterion(icond) = crit1;
        
        %         end
    end
    %     optimal_criterion = mean(optimal_criterion,1);
    saved_pars = pars;
end
dec_rule = optimal_criterion;


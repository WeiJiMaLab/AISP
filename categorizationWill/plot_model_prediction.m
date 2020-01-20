function plot_model_prediction(type)

data = getData();
switch type
    case 'Bayes'
        p = load('pars/parsBayes.mat');
        l = p.likelihoods;
        [~,maxID] = min(l,[],2);
        pars = zeros(size(p.pars,1),size(p.pars,2));
        for iSubj = 1:size(p.pars,1)
            pars(iSubj,:) = p.pars(iSubj,:,maxID(iSubj));
            sigmas = exp(pars(iSubj,1:6));
            beta = [pars(iSubj,7),exp(pars(iSubj,8))];
            lambda = pars(iSubj,9);
            data(data(:,1)==iSubj,6) = 2*bayes_simulate(data(data(:,1)==iSubj,2:3),sigmas,beta,lambda)-1;
        end
    case 'Freq'
        p = load('pars/parsFreq.mat');
        l = p.likelihoods;
        [~,maxID] = min(l,[],2);
        pars = zeros(size(p.pars,1),size(p.pars,2));
        for iSubj = 1:size(p.pars,1)
            pars(iSubj,:) = p.pars(iSubj,:,maxID(iSubj));
            sigmas = exp(pars(iSubj,1:6));
            beta = [pars(iSubj,7),exp(pars(iSubj,8))];
            lambda = pars(iSubj,9);
            data(data(:,1)==iSubj,6) = 2*freq_simulate(data(data(:,1)==iSubj,2:3),sigmas,beta,lambda)-1;
        end
    case 'Freq2'
        p = load('pars/parsFreq2.mat');
        l = p.likelihoods;
        [~,maxID] = min(l,[],2);
        pars = zeros(size(p.pars,1),size(p.pars,2));
        for iSubj = 1:size(p.pars,1)
            pars(iSubj,:) = p.pars(iSubj,:,maxID(iSubj));
            sigmas = exp(pars(iSubj,1:6));
            beta = [pars(iSubj,7),exp(pars(iSubj,8))];
            lambda = pars(iSubj,9);
            data(data(:,1)==iSubj,6) = 2*freq2_simulate(data(data(:,1)==iSubj,2:3),sigmas,beta,lambda)-1;
        end
    case 'Freq3'
        p = load('pars/parsFreq3.mat');
        l = p.likelihoods;
        [~,maxID] = min(l,[],2);
        pars = zeros(size(p.pars,1),size(p.pars,2));
        for iSubj = 1:size(p.pars,1)
            pars(iSubj,:) = p.pars(iSubj,:,maxID(iSubj));
            sigmas = exp(pars(iSubj,1:6));
            beta = [pars(iSubj,7),exp(pars(iSubj,8))];
            lambda = pars(iSubj,9);
            data(data(:,1)==iSubj,6) = 2*freq3_simulate(data(data(:,1)==iSubj,2:3),sigmas,beta,lambda)-1;
        end
    case 'Var'
        p = load('pars/parsVar.mat');
        l = p.likelihoods;
        [~,maxID] = min(l,[],2);
        pars = zeros(size(p.pars,1),size(p.pars,2));
        for iSubj = 1:size(p.pars,1)
            pars(iSubj,:) = p.pars(iSubj,:,maxID(iSubj));
            sigmas = exp(pars(iSubj,1:6));
            beta = [pars(iSubj,7),exp(pars(iSubj,8))];
            lambda = pars(iSubj,9);
            data(data(:,1)==iSubj,6) = 2*var_simulate(data(data(:,1)==iSubj,2:3),sigmas,beta,lambda)-1;
        end
end


Nbins = 11;
Nbins2 = 100;

for iSubj = 1:11
    datSubj = data(data(:,1)==iSubj,:);
    figure
    subplot(1,2,1)
    bins = quantile(datSubj(:,3),linspace(0,1,Nbins+1));
    x = zeros(Nbins,6);
    y = zeros(Nbins,6);
    yModel = zeros(Nbins,6);
    for k = 1:Nbins
        d = datSubj((datSubj(:,3)<bins(k+1)) & (datSubj(:,3)>bins(k)),:);
        for iNoise = 1:6
            dNoise = d(d(:,2)==iNoise,:);
            x(k,iNoise) = mean(dNoise(:,3));
            y(k,iNoise) = 0.5+0.5*mean(dNoise(:,4));
            yModel(k,iNoise) = 0.5+0.5*mean(dNoise(:,6));
        end
    end
    plot(x,y,'k.-','MarkerSize',10)
    hold on
    plot(x,yModel,'bs--','MarkerSize',10)
    xlabel('Stimulus Level','FontSize',16)
    ylabel('Proportion Category 1','FontSize',16)
    ylim([0,1])
    box off
    set(gca,'TickDir','out')
    
    subplot(1,2,2)
    k = 0;
    x = zeros(Nbins2,1);
    y = zeros(Nbins2,1);
    yModel = zeros(Nbins2,1);
    for i = linspace(-40,40,Nbins2)
        k = k+1;
        d = datSubj((datSubj(:,3)<i) & (datSubj(:,3)>(i-1)),:);
        x(k) = i-0.5;
        y(k) = 0.5+0.5*mean(d(:,4));
        yModel(k) = 0.5+0.5*mean(d(:,6));
    end
    plot(x,y,'k.','MarkerSize',10)
    hold on
    plot(x,yModel,'b.','MarkerSize',10)
    xlabel('Stimulus Level','FontSize',16)
    ylabel('Proportion Category 1','FontSize',16)
    ylim([0,1])
    box off
    set(gca,'TickDir','out')
    set(gcf,'Position',[560,723,560,225])
end



datSubj = data;
figure
subplot(1,2,1)
bins = quantile(datSubj(:,3),linspace(0,1,Nbins+1));
x = zeros(Nbins,6);
y = zeros(Nbins,6);
yModel = zeros(Nbins,6);
for k = 1:Nbins
    d = datSubj((datSubj(:,3)<bins(k+1)) & (datSubj(:,3)>bins(k)),:);
    for iNoise = 1:6
        dNoise = d(d(:,2)==iNoise,:);
        x(k,iNoise) = mean(dNoise(:,3));
        y(k,iNoise) = 0.5+0.5*mean(dNoise(:,4));
        yModel(k,iNoise) = 0.5+0.5*mean(dNoise(:,6));
    end
end
plot(x,y,'k.-','MarkerSize',10)
hold on
plot(x,yModel,'bs--','MarkerSize',10)
xlabel('Stimulus Level','FontSize',16)
ylabel('Proportion Category 1','FontSize',16)
ylim([0,1])
box off
set(gca,'TickDir','out')

subplot(1,2,2)
k = 0;
x = zeros(Nbins2,1);
y = zeros(Nbins2,1);
yModel = zeros(Nbins2,1);
for i = linspace(-40,40,Nbins2)
    k = k+1;
    d = datSubj((datSubj(:,3)<i) & (datSubj(:,3)>(i-1)),:);
    x(k) = i-0.5;
    y(k) = 0.5+0.5*mean(d(:,4));
    yModel(k) = 0.5+0.5*mean(d(:,6));
end
plot(x,y,'k.','MarkerSize',10)
hold on
plot(x,yModel,'b.','MarkerSize',10)
xlabel('Stimulus Level','FontSize',16)
ylabel('Proportion Category 1','FontSize',16)
ylim([0,1])
box off
set(gca,'TickDir','out')
set(gcf,'Position',[560,723,560,225])

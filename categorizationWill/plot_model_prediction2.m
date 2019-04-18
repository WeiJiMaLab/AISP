function plot_model_prediction2(type)

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
end


for iSubj = 1:11
    d = data(data(:,1)==iSubj,:);
    figure
    x = 1:6;
    y = zeros(2,6);
    yModel = zeros(2,6);
    for iNoise = 1:6
        dNoise = d(d(:,2)==iNoise,:);
        y(1,iNoise) = 0.5+0.5*mean(dNoise(dNoise(:,5)==1,4));
        y(2,iNoise) = 0.5+0.5*mean(dNoise(dNoise(:,5)==-1,4));
        yModel(1,iNoise) = 0.5+0.5*mean(dNoise(dNoise(:,5)==1,6));
        yModel(2,iNoise) = 0.5+0.5*mean(dNoise(dNoise(:,5)==-1,6));
    end
    plot(x,y(1,:),'b.-','MarkerSize',10)
    hold on
    plot(x,y(2,:),'r.-','MarkerSize',10)
    plot(x,yModel(1,:),'bs--','MarkerSize',10)
    plot(x,yModel(2,:),'rs--','MarkerSize',10)
    xlabel('Reliability Level','FontSize',16)
    ylabel('Proportion Category 1','FontSize',16)
    ylim([0,1])
    xlim([0,7])
    box off
    set(gca,'TickDir','out')
end



d = data;
figure
x = 1:6;
y = zeros(2,6);
yModel = zeros(2,6);
for iNoise = 1:6
    dNoise = d(d(:,2)==iNoise,:);
    y(1,iNoise) = 0.5+0.5*mean(dNoise(dNoise(:,5)==1,4));
    y(2,iNoise) = 0.5+0.5*mean(dNoise(dNoise(:,5)==-1,4));
    yModel(1,iNoise) = 0.5+0.5*mean(dNoise(dNoise(:,5)==1,6));
    yModel(2,iNoise) = 0.5+0.5*mean(dNoise(dNoise(:,5)==-1,6));
end
plot(x,y(1,:),'b.-','MarkerSize',10)
hold on
plot(x,y(2,:),'r.-','MarkerSize',10)
plot(x,yModel(1,:),'bs--','MarkerSize',10)
plot(x,yModel(2,:),'rs--','MarkerSize',10)
xlabel('Reliability Level','FontSize',16)
ylabel('Proportion Category 1','FontSize',16)
ylim([0,1])
xlim([0,7])
box off
set(gca,'TickDir','out')
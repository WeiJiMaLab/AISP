data = getData();

Nbins = 11;
Nbins2 = 100;

for iSubj = 1:11
    datSubj = data(data(:,1)==iSubj,:);
    figure
    subplot(1,2,1)
    bins = quantile(datSubj(:,3),linspace(0,1,Nbins+1));
    x = zeros(Nbins,6);
    y = zeros(Nbins,6);
    for k = 1:Nbins
        d = datSubj((datSubj(:,3)<bins(k+1)) & (datSubj(:,3)>bins(k)),:);
        for iNoise = 1:6
            dNoise = d(d(:,2)==iNoise,:);
            x(k,iNoise) = mean(dNoise(:,3));
            y(k,iNoise) = 0.5+0.5*mean(dNoise(:,4));
        end
    end
    plot(x,y,'.-','MarkerSize',10)
    hold on
    xlabel('Stimulus Level','FontSize',16)
    ylabel('Proportion Category 1','FontSize',16)
    ylim([0,1])
    box off
    set(gca,'TickDir','out')
    
    subplot(1,2,2)
    k = 0;
    x = zeros(Nbins2,1);
    y = zeros(Nbins2,1);
    for i = linspace(-40,40,Nbins2)
        k = k+1;
        d = datSubj((datSubj(:,3)<i) & (datSubj(:,3)>(i-1)),:);
        x(k) = i-0.5;
        y(k) = 0.5+0.5*mean(d(:,4));
    end
    plot(x,y,'k.','MarkerSize',10)
    hold on
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
for k = 1:Nbins
    d = datSubj((datSubj(:,3)<bins(k+1)) & (datSubj(:,3)>bins(k)),:);
    for iNoise = 1:6
        dNoise = d(d(:,2)==iNoise,:);
        x(k,iNoise) = mean(dNoise(:,3));
        y(k,iNoise) = 0.5+0.5*mean(dNoise(:,4));
    end
end
plot(x,y,'.-','MarkerSize',10)
hold on
xlabel('Stimulus Level','FontSize',16)
ylabel('Proportion Category 1','FontSize',16)
ylim([0,1])
box off
set(gca,'TickDir','out')

subplot(1,2,2)
k = 0;
x = zeros(Nbins2,1);
y = zeros(Nbins2,1);
for i = linspace(-40,40,Nbins2)
    k = k+1;
    d = datSubj((datSubj(:,3)<i) & (datSubj(:,3)>(i-1)),:);
    x(k) = i-0.5;
    y(k) = 0.5+0.5*mean(d(:,4));
end
plot(x,y,'k.','MarkerSize',10)
hold on
xlabel('Stimulus Level','FontSize',16)
ylabel('Proportion Category 1','FontSize',16)
ylim([0,1])
box off
set(gca,'TickDir','out')
set(gcf,'Position',[560,723,560,225])

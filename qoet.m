%% Simulate
% function [] = qoeareal()
clear;
close all force;
% profile on
% lambda = 10;
% mu = 15;
% mu = 10;
% lambda = [0.5 * mu, 0.8*mu];

arrival = 'M';
service = 'M';
% arrival and service
% M - Markovian
% G - Gerneral distribution
% AR - auto correlation with lag 1 and 0.8
% policy = [1 2 3];
policy = 2;
% policies:
% 1 - n-policy
% 2 - D-policy
% 3 - T-policy
v=2;
a = [0.8 1 1.2];
mu = 1;
lambda = a .* mu;
simtime = [30 60 120 240 480 960 1920 3840];
% a=0.5;
n = length(simtime);
o = length(a);
m = 100; % runs
QoE = zeros(m,o,n);
E1vod = zeros(m,o,n);
E2vod = zeros(m,o,n);
E1live = zeros(m,o,n);
E2live = zeros(m,o,n);
L = zeros(m,o,n);
N = zeros(m,o,n);
rate = zeros(m,o,n);
ratio = zeros(m,o,n);
StallingRatio = zeros(m,o,n);
reallambda = zeros(m,o,n);
realmu = zeros(m,o,n);
reala = zeros(m,o,n);
realcv1 = zeros(m,o,n);
realcv2 = zeros(m,o,n);
permutevideo = 1;
circshiftvideo = 0;
load framesizeData;
% load timeOfFrameDownload;
% Sij = frameSize(:,3);
Sij = sum(frameSize'); % video bitrate / 24
% Sij=Sij(end:-1:1);
    dmax = 5; % sec
    d = dmax * mu; % frames
tic
waitb = waitbar(0,'1','Name','Simulating Video Streaming');
for j=1:m
%     [bwPerSecond1,Sijx]=getData(j-1,1);
%     lambdaeingabe = bwPerSecond1 / 8 * 1000; % bandwidth
%     lambda=mean(lambdaeingabe);
%     if permutevideo == 1
%         idx = randperm(length(Sij)); %shuffle framess
%         Sij = Sij(idx);
%     end
%     if circshiftvideo == 1
%         Sij = circshift(Sij,[0 ceil(length(Sij)/7)]);
%     end
% %     Sij = Sij(:,3);
%     mu = mean(Sij)*24;
%     a0 = lambda./mu;
%     
    

%     lambdaeingabe = 
    
    % dmax = 1:0.1:20;
    % lambda = 1:24;
    % v = [1 2 3]; % variations coefficent
    for k=1:o
%         timeOfFrameDownload = 0;
%         timeOfFrameDownload = getTimeOfFrameDownload(lambdaeingabe/a0*a(k), Sij);
%         timeOfFrameDownload{j,k} = sim4(lambdaeingabe, Sij(1:end)*a0/a(k));
        for l=1:n
%             tic
            waitbar(j/m,waitb,['traffic pattern ' num2str(j) '/' num2str(m) ', reception ratio ' num2str(k) '/' num2str(o) ', policy ' num2str(l) '/' num2str(n)]);
%             [QoE(j,k,l),E1vod(j,k,l),E2vod(j,k,l),E1live(j,k,l),E2live(j,k,l),N(j,k,l),L(j,k,l),rate(j,k,l),ratio(j,k,l),StallingRatio(j,k,l)] = runSim3withData(lambdaeingabe/a0*a(k),Sij,timeOfFrameDownload,dmax,v(l),0,arrival,service,policy);
            [QoE(j,k,l),E1vod(j,k,l),E2vod(j,k,l),E1live(j,k,l),E2live(j,k,l),N(j,k,l),L(j,k,l),rate(j,k,l),ratio(j,k,l),StallingRatio(j,k,l),reallambda(j,k,l), realmu(j,k,l), reala(j,k,l),realcv1(j,k,l),realcv2(j,k,l)] = runSim2(lambda(k),mu,dmax,v,0,arrival,service,policy,simtime(l));
%             toc
        end
    end
    % v=repmat(v',1,3);
end
toc
close(waitb);

AnaL = dmax.*mu./lambda;
AnaN = (1-lambda./mu)./dmax;
AnaQoE = exp(-(0.15 .* AnaL + 0.19) .* AnaN);
AnaQoEplus = exp(-0.15 .* AnaL .* AnaN) + exp(- 0.19 .* AnaN) - 1;


[nL, nN, nQoE] = npolicy(mu,lambda,d);
AnaQoE = (AnaQoE * 3.5) +1.5;
QoE = (QoE * 3.5) + 1.5;
for l=1:o
    x{l} = ['\lambda = ' num2str(a(l)) ' Mbit/s'];
end
% x={'n-policy','D-policy','T-policy'};
save(['results/simtime' arrival service '_3']);
%% Plot results
close all;
conffactor = 1.96;
% clear x;
for l=1:o
    x{l} = ['\lambda = ' num2str(a(l)) ' Mbit/s'];
end
% x={'n-policy','D-policy','T-policy'};

figure(2)
box on;
% mx = reshape(mean(QoE,2),n,o);
% stdx = reshape(std(QoE,0,2),n,o);
% for k=1:o
%     errorbar(v,mx(:,k),stdx(:,k));
%     hold on;
% end
mQoE = reshape(nanmean(QoE),o,n);
sQoE = conffactor * reshape(nanstd(QoE),o,n) / sqrt(m);
errorbar(simtime,mQoE(1,:),sQoE(1,:),'-k','LineWidth',2);
hold on;
errorbar(simtime,mQoE(2,:),sQoE(2,:),'--k','LineWidth',2);
errorbar(simtime,mQoE(3,:),sQoE(3,:),':k','LineWidth',2);
set(gca,'XScale','log');
% xlim([v(1) v(end)])
ylim([1.5 5])
% plot(v,QoE);
hold on;
% plot(v,AnaQoE,'k')
% plot([v(1) v(end)],[AnaQoE; AnaQoE],'k')
hold on;
ylabel('mean opinion score MOS')
xlabel('coefficient of variation c_v')
legend(x,'Location','SouthEast');
set(gca,'Fontsize',14)
saveas(gcf,['figs\corr' arrival service '_QoE_' int2str(permutevideo) int2str(circshiftvideo)],'eps2c');

figure(4)
box on
set(gca,'Fontsize',14)
mN = reshape(nanmean(N),o,n);
sN = conffactor * reshape(nanstd(N),o,n) / sqrt(m);
errorbar(v,mN(1,:),sN(1,:),'-k','LineWidth',2);
hold on;
errorbar(v,mN(2,:),sN(2,:),'--k','LineWidth',2);
errorbar(v,mN(3,:),sN(3,:),':k','LineWidth',2);
xlim([v(1) v(end)])
ylim([0 max(get(gca,'YLim'))])
% plot(v,QoE);
hold on;
% plot(v,AnaQoE,'k')
% plot([v(1) v(end)],[AnaQoE; AnaQoE],'k')
hold on;
ylabel('mean frequency of stalling events N')
xlabel('coefficient of variation c_v')
legend(x,'Location','NorthEast');
set(gca,'Fontsize',14)
saveas(gcf,['figs\corr' arrival service '_N_' int2str(permutevideo) int2str(circshiftvideo)],'eps2c');

L(L==0)=NaN;
figure(5)
box on
set(gca,'Fontsize',14)
mL = reshape(nanmean(L),o,n);
sL = conffactor * reshape(nanstd(L),o,n) / sqrt(m);
errorbar(v,mL(1,:),sL(1,:),'-k','LineWidth',2);
hold on;
errorbar(v,mL(2,:),sL(2,:),'--k','LineWidth',2);
errorbar(v,mL(3,:),sL(3,:),':k','LineWidth',2);
xlim([v(1) v(end)])
ylim([0 100])
% plot(v,QoE);
hold on;
% ylim([0 20]);
% plot(v,AnaQoE,'k')
% plot([v(1) v(end)],[AnaQoE; AnaQoE],'k')
hold on;
ylabel('mean length of stalling events L')
xlabel('coefficient of variation c_v')
legend(x,'Location','NorthEast');
set(gca,'Fontsize',14)
saveas(gcf,['figs\corr' arrival service '_L_' int2str(permutevideo) int2str(circshiftvideo)],'eps2c');

figure(61)
cdfplot(reallambda(:))
xlabel('mean bandwidth of a run in Mbit/s')

figure(62)
cdfplot(realmu(:))
xlabel('mean video bit rate of a run in Mbit/s')

figure(63)
cdfplot(realcv1(:))
xlabel('c_v of the bandwidth')

figure(64)
cdfplot(realcv2(:))
xlabel('c_v of the bit rate')

% figure(22)
% errorbar(v,mean(QoE'),std(QoE'));
% xlim([0 v(end)])
% ylim([0 1])
% % plot(v,QoE);
% hold on;
% % plot(v,AnaQoE,'k')
% plot([v(1) v(end)],[AnaQoEplus; AnaQoEplus],'k')
% hold on;
% ylabel('QoE value')
% xlabel('v')
% 
% figure(4)
% 
% errorbar(repmat(a',1,n),reshape(mean(N),o,n),conffactor * reshape(std(N),o,n) / sqrt(m));
% xlim([0 a(end)])
% ylim([0 max(max(max(N)))])
% % plot(v,N)
% hold on;
% % plot(v,AnaN,'k')
% % plot([v(1) v(end)],[AnaN; AnaN],'k')
% hold on;
% ylabel('N^*')
% xlabel('a')
% % legend(x);
% saveas(gcf,['figs\areal' arrival service num2str(policy) '_N'],'eps2c');
% 
% figure(6)
% 
% errorbar(repmat(a',1,n),reshape(mean(StallingRatio),o,n),conffactor * reshape(std(StallingRatio),o,n) / sqrt(m));
% xlim([0 a(end)])
% ylim([0 max(max(max(StallingRatio)))])
% % plot(v,N)
% hold on;
% % plot(v,AnaN,'k')
% % plot([v(1) v(end)],[1-a; 1-a],'k')
% hold on;
% ylabel('stalling ratio: ')
% % totalstalling / (totalplaying + totalstalling)
% xlabel('a')
% % legend(x);
% saveas(gcf,['figs\areal' arrival service num2str(policy) '_ratio'],'eps2c');
% 
% figure(5)
% 
% errorbar(repmat(a',1,n),reshape(mean(L),o,n),conffactor * reshape(std(L),o,n) / sqrt(m));
% xlim([0 a(end)])
% ylim([0 max(max(max(L)))])
% % semilogy(v,L)
% hold on;
% % semilogy(v,AnaL,'k')
% % plot([v(1) v(end)],[AnaL; AnaL],'k')
% hold on;
% ylabel('total stalling duration')
% xlabel('a')
% % legend(x);
% saveas(gcf,['figs\areal' arrival service num2str(policy) '_L'],'eps2c');
% %
% % load('fittings');
% % x1 = pLIVE(1) * exp(-pLIVE(2) * rate) + pLIVE(3) * exp(-pLIVE(4) * rate);
% % x2 = pVOD(1) * exp(-pVOD(2) * rate) + pVOD(3) * exp(-pVOD(4) * rate);
% % x3 = qLIVE(1) * exp(-qLIVE(2) * StallingRatio) + qLIVE(3) * exp(-qLIVE(4) * StallingRatio);
% % x4 = qVOD(1) * exp(-qVOD(2) * StallingRatio) + qVOD(3) * exp(-qVOD(4) * StallingRatio);
% %
% % figure(71)
% % errorbar(v,mean(x1'),std(x1'));
% % hold on;
% % errorbar(v,mean(x3'),std(x3'));
% % xlabel('v')
% % ylabel('Play time (min)')
% % legend({'ratio','rate'})
% % % figure(72)
% % errorbar(v,mean(x2'),std(x2'));
% % hold on
% % errorbar(v,mean(x4'),std(x4'));
% % xlabel('v')
% % ylabel('Play time (min)')
% % legend({'ratio','rate'})
% %
% % figure(3)
% % errorbar(v,reshape(mean(E1vod,2),n,o),reshape(std(E1vod,0,2),n,o));
% % hold all;
% % errorbar(v,reshape(mean(E2vod,2),n,o),reshape(std(E2vod,0,2),n,o));
% % errorbar(v,reshape(mean(E1live,2),n,o),reshape(std(E1live,0,2),n,o));
% % errorbar(v,reshape(mean(E2live,2),n,o),reshape(std(E2live,0,2),n,o));
% % ylabel('Play time [min]');
% % xlabel('v')
% % legend({'rateVOD','ratioVOD','rateLIVE','ratioLIVE'},'Location','SouthEast')
% % saveas(gcf,[arrival service num2str(policy) '_fit'],'eps2c');
% % xlim([0 3.2]);
%%
% only a = 0.5
figure(3)
box on
set(gca,'Fontsize',14)
mE1vod = reshape(mean(E1vod),o,n);
sE1vod = conffactor * reshape(std(E1vod),o,n) / sqrt(m);
% errorbar(a,mE1vod(:,1),sE1vod(:,1),'-k','LineWidth',2);
hold all;
% errorbar(a,mE1vod(:,2),sE1vod(:,2),'--k','LineWidth',2);
% errorbar(a,mE1vod(:,3),sE1vod(:,3),':k','LineWidth',2);
mE2vod = reshape(mean(E2vod),o,n);
sE2vod = conffactor * reshape(std(E2vod),o,n) / sqrt(m);
errorbar(a,mE2vod(:,1),sE2vod(:,1),'-k','LineWidth',2);
errorbar(a,mE2vod(:,2),sE2vod(:,2),'--k','LineWidth',2);
errorbar(a,mE2vod(:,3),sE2vod(:,3),':k','LineWidth',2);
xlim([min(a),max(a)])
ylim([0 30])
% errorbar(a,mean(E1live),conffactor * std(E1live) / sqrt(m),'-r');
% errorbar(a,mean(E2live),conffactor * std(E2live) / sqrt(m),'--r');
ylabel('play time [min]');
xlabel('offered load a')
legend(x,'Location','SouthEast')
% legend({'rateVOD','ratioVOD'})
saveas(gcf,['figs\corr' arrival service '_fit_' int2str(permutevideo) int2str(circshiftvideo)],'eps2c');

% figure(6)
% plot(lambda/mu,rate)
% ylabel('rate');
%
% figure(7)
% plot(lambda/mu,ratio)
% ylabel('ratio');

% figure(4)
% [h1,stats1] = cdfplot(QoE(:,11,1));
% hold all;
% [h2,stats2] = cdfplot(QoE(:,11,2));
% [h3,stats3] = cdfplot(QoE(:,11,3));
% legend(x,'Location','NorthWest')
% xlim([1.5 5])
% h1.LineStyle = '-';
% h2.LineStyle = '--';
% h3.LineStyle = ':';
% h1.LineWidth = 2;
% h2.LineWidth = 2;
% h3.LineWidth = 2;
% % h1.Color = 'k';
% % h2.Color = 'k';
% % h3.Color = 'k';
% xlabel('MOS at a = 1')
% saveas(gcf,['figs\areal' arrival service '_a_' int2str(permutevideo) int2str(circshiftvideo)],'eps2c');
%%
figure(5)
scatter(1:30,QoE(:,11,1));
hold all;
scatter(1:30,QoE(:,11,2),'s');
scatter(1:30,QoE(:,11,3),'*');
legend(x,'Location','NorthWest')
ylabel('MOS at a = 1')
xlabel('run with shift of traffic')
saveas(gcf,['figs\areal' arrival service '_scatter_' int2str(permutevideo) int2str(circshiftvideo)],'eps2c');

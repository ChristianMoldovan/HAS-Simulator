
load results/MM_3
close all
conffactor = 1.96;
% figure(1)

% Number of runs with stalling
figure(2)
box on
set(gca,'Fontsize',14)
N1 = squeeze(sum(N>0));
% mN = reshape(nanmean(N),o,n);
% sN = conffactor * reshape(nanstd(N),o,n) / sqrt(m);
% errorbar(v,mN(1,:),sN(1,:),'-k','LineWidth',2);
% hold on;
% errorbar(v,mN(2,:),sN(2,:),'--k','LineWidth',2);
% errorbar(v,mN(3,:),sN(3,:),':k','LineWidth',2);
cols = copper(3);
plot(v,N1(1,:),'-k','LineWidth',2,'Color',cols(1,:));
hold on;
plot(v,N1(2,:),'--k','LineWidth',2,'Color',cols(2,:));
plot(v,N1(3,:),':k','LineWidth',2,'Color',cols(3,:));
xlim([v(1) v(end)])
ylim([0 max(size(N))+3])
% plot(v,QoE);
hold on;
% plot(v,AnaQoE,'k')
% plot([v(1) v(end)],[AnaQoE; AnaQoE],'k')
hold on;
ylabel('Number of runs where stalling occured')
xlabel('coefficient of variation c_v')
legend(x,'Location','NorthEast');
set(gca,'Fontsize',14)
set(gca,'XScale','log')
saveas(gcf,['figs\r1' arrival service '_N_' int2str(permutevideo) int2str(circshiftvideo)],'eps2c');


% QoE vs measured a
figure(3)
plot(squeeze(reala(:,1,:)),squeeze(QoE(:,1,:)),'ko')
hold on;
plot(squeeze(reala(:,2,:)),squeeze(QoE(:,2,:)),'r.')
plot(squeeze(reala(:,3,:)),squeeze(QoE(:,3,:)),'bs')
xlim([0.5 1.5])
ylim([1.5 5])
% plot(v,QoE);
hold on;
% plot(v,AnaQoE,'k')
% plot([v(1) v(end)],[AnaQoE; AnaQoE],'k')
hold on;
ylabel('mean opinion score MOS')
xlabel('measured bandwidth in Mbit/s')
% legend(x,'Location','SouthEast');
set(gca,'Fontsize',14)
saveas(gcf,['figs\r1' arrival service '_QoE_' int2str(permutevideo) int2str(circshiftvideo)],'eps2c');

%%

load results/simtimeMG_3

close all;
conffactor = 1.96;
% clear x;
% for l=1:o
%     x{l} = ['\lambda = ' num2str(a(l)) ' Mbit/s'];
% end
% x={'n-policy','D-policy','T-policy'};
fontsize = 24;
figure(2)
box on;
% mx = reshape(mean(QoE,2),n,o);
% stdx = reshape(std(QoE,0,2),n,o);
% for k=1:o
%     errorbar(v,mx(:,k),stdx(:,k));
%     hold on;
% end
cols = copper(3);
mQoE = reshape(nanmean(QoE),o,n);
sQoE = conffactor * reshape(nanstd(QoE),o,n) / sqrt(m);
errorbar(simtime,mQoE(1,:),sQoE(1,:),'-','LineWidth',2,'Color',cols(1,:));
hold on;
errorbar(simtime,mQoE(2,:),sQoE(2,:),'--','LineWidth',2,'Color',cols(2,:));
errorbar(simtime,mQoE(3,:),sQoE(3,:),':','LineWidth',2,'Color',cols(3,:));
set(gca,'XScale','log');
xlim([simtime(1)*0.9 simtime(end)*1.1])
ylim([1.5 5])
% plot(v,QoE);
hold on;
% plot(v,AnaQoE,'k')
% plot([v(1) v(end)],[AnaQoE; AnaQoE],'k')
hold on;
ylabel('mean opinion score MOS')
xlabel('simulation time [s]')
legend(x,'Location','East');
set(gca,'Fontsize',fontsize)
saveas(gcf,['figs\simtime_' arrival service '_QoE_' int2str(permutevideo) int2str(circshiftvideo)],'eps2c');

figure(3)
box on;
% mx = reshape(mean(QoE,2),n,o);
% stdx = reshape(std(QoE,0,2),n,o);
% for k=1:o
%     errorbar(v,mx(:,k),stdx(:,k));
%     hold on;
% end
% mN = reshape(nanmean(N),o,n);
% sN = conffactor * reshape(nanstd(N),o,n) / sqrt(m);
N1 = squeeze(sum(N>0))/100;
sN1 = 1.96*sqrt(N1.*(1-N1)/100);
errorbar(simtime,N1(1,:),sN1(1,:),'-','LineWidth',2,'Color',cols(1,:));
hold on;
errorbar(simtime,N1(2,:),sN1(2,:),'--','LineWidth',2,'Color',cols(2,:));
errorbar(simtime,N1(3,:),sN1(3,:),':','LineWidth',2,'Color',cols(3,:));
set(gca,'XScale','log');
% xlim([simtime(1)*0.9 simtime(end)*1.1])
% ylim([1.5 5])
% plot(v,QoE);
hold on;
% plot(v,AnaQoE,'k')
% plot([v(1) v(end)],[AnaQoE; AnaQoE],'k')
hold on;
xlim([simtime(1)*0.9 simtime(end)*1.1])
ylim([0 1])
ylabel('ratio of runs with stalling')
xlabel('simulation time [s]')
legend(x,'Location','SouthEast');
if (strcmp(arrival,'M'))
    legend(x,'Location','East');
end
set(gca,'Fontsize',fontsize)
saveas(gcf,['figs\simtime_' arrival service '_N_' int2str(permutevideo) int2str(circshiftvideo)],'eps2c');


%% Simulate
% function [] = qoeareal()
clear;
close all force;
% profile on
% lambda = 10;
% mu = 15;
% mu = 10;
% lambda = [0.5 * mu, 0.8*mu];

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
% v=logspace(-1,1,10);
arrival = 'MG';
service = 'MG';
v = [0.1 1 10];
a = [0.8 1 1.2];
simtime = [30 150 750];
mu = 1;
lambda = a .* mu;
% a=0.5;
q = length(arrival);
r = length(service);
n = length(v);
o = length(a);
p = length(simtime);
m = 100; % runs
QoE = zeros(m,o,n,p,q,r);
E1vod = zeros(m,o,n,p,q,r);
E2vod = zeros(m,o,n,p,q,r);
E1live = zeros(m,o,n,p,q,r);
E2live = zeros(m,o,n,p,q,r);
L = zeros(m,o,n,p,q,r);
N = zeros(m,o,n,p,q,r);
rate = zeros(m,o,n,p,q,r);
ratio = zeros(m,o,n,p,q,r);
StallingRatio = zeros(m,o,n,p,q,r);
reallambda = zeros(m,o,n,p,q,r);
realmu = zeros(m,o,n,p,q,r);
reala = zeros(m,o,n,p,q,r);
realcv1 = zeros(m,o,n,p,q,r);
realcv2 = zeros(m,o,n,p,q,r);
permutevideo = 1;
circshiftvideo = 0;
% load framesizeData;
% load timeOfFrameDownload;
% Sij = frameSize(:,3);
% Sij = sum(frameSize'); % video bitrate / 24
% Sij=Sij(end:-1:1);
dmax = 5; % sec
d = dmax * mu; % frames
tic
waitb = waitbar(0,'1','Name','Simulating Video Streaming');

for j=1:m % 100 runs
    for k=1:o % network bandwidth
        for l=1:n % network variation
            for h=1:p         % simtime
                for x=1:q % arrival
                    for y=1:r % service
                        %             tic
                        waitbar(j/m,waitb,['traffic pattern ' num2str(j) '/' num2str(m) ', reception ratio ' num2str(k) '/' num2str(o) ', policy ' num2str(l) '/' num2str(n)]);
                        [QoE(j,k,l,h,x,y),E1vod(j,k,l,h,x,y),E2vod(j,k,l,h,x,y),E1live(j,k,l,h,x,y),E2live(j,k,l,h,x,y),N(j,k,l,h,x,y),L(j,k,l,h,x,y),rate(j,k,l,h,x,y),ratio(j,k,l,h,x,y),StallingRatio(j,k,l,h,x,y),reallambda(j,k,l,h,x,y), realmu(j,k,l,h,x,y), reala(j,k,l,h,x,y),realcv1(j,k,l,h,x,y),realcv2(j,k,l,h,x,y)] = runSim2(lambda(k),mu,dmax,v(l),0,arrival(x),service(y),policy,simtime(h));
                        %             toc
                    end
                end
            end
        end
    end
end
toc
close(waitb);

% AnaL = dmax.*mu./lambda;
% AnaN = (1-lambda./mu)./dmax;
% AnaQoE = exp(-(0.15 .* AnaL + 0.19) .* AnaN);
% AnaQoEplus = exp(-0.15 .* AnaL .* AnaN) + exp(- 0.19 .* AnaN) - 1;
% 
% 
% [nL, nN, nQoE] = npolicy(mu,lambda,d);
% AnaQoE = (AnaQoE * 3.5) +1.5;
QoE = (QoE * 3.5) + 1.5;
% for l=1:o
%     x{l} = ['\lambda = ' num2str(a(l)) ' Mbit/s'];
% end
% x={'n-policy','D-policy','T-policy'};
save(['results/doe_' int2str(m)]);
%%
close all
figure(1)
% set(gca,'FontSize',24)
qoesize = size(QoE);
repa = repmat(a,qoesize(1),1,qoesize(3),qoesize(4),qoesize(5),qoesize(6));

repsimtime = zeros(qoesize);
repsimtime(:,:,:,1,:,:)=simtime(1);
repsimtime(:,:,:,2,:,:)=simtime(2);
repsimtime(:,:,:,3,:,:)=simtime(3);

repv = zeros(qoesize);
repv(:,:,1,:,:,:)=v(1);
repv(:,:,2,:,:,:)=v(2);
repv(:,:,3,:,:,:)=v(3);

% [w,b]=histc(simtimeagg,[0 40 150 600 inf]);
[w1,realcv1x]=histc(realcv1,[0.4 2 inf]);
[w2,realcv2x]=histc(realcv2,[0.4 2 inf]);
% [w,h]=histc(realaagg,[0.8 0.9 1 1.1 1.2 inf]);
fontsize = 12;
% figure(1)
% set(gca,'Fontsize',fontsize);
h=maineffectsplot(QoE(:),{repa(:),repsimtime(:),realcv1x(:),realcv2x(:)},'varnames',{'\lambda','simtime','measured c_\lambda','measured c_\mu'});
for i=1:length(h.Children)
    set(h.Children(i),'FontSize',fontsize);
    h.Children(i).Children(1).LineWidth = 3;
end
% set(h.Children(5).YLabel.String,'MOS');
h.Children(end).YLabel.String = 'MOS';

saveas(gcf,['figs\doe_' int2str(m)],'eps2c');




function [QoE,E1vod,E2vod,E1live,E2live,N,L,rate,ratio,StallingRatio] = runSim3withData(bwPerSecond,Sij,timeOfFrameDownload,d,cv,plotme,arrival,service,policy,ar)
% close all;
% lambda = 10; % fps
% mu = 25; % fps
% close all;
% if nargin<6
%     plotme = false;
% end
dmin = 0;
% dmax = d;
% dmax = 50; % f
% Sij contains 2 sec

if nargin <= 9
    ar = 0.89;
end

lambda = mean(bwPerSecond);
% lambda = 1/mean(diff(timeOfFrameDownload));
mu = mean(Sij);
a = lambda/mu/24;
duration = 100000;
fps = 24;
framesize = 1/fps; % s
if (strcmp(service,'D'))
    tmu = Sij;
else
    tmu = generateVideo(service, duration, cv, mu,ar);
end
Sij = tmu;

% if (plotme==1)
%     timeOfFrameDownload = getTimeOfFrameDownload(bwPerSecond, Sij);
% end

if (strcmp(arrival,'D'))
%     tlambda = bwPerSecond;
    tlambda = diff(timeOfFrameDownload);
% tlambda = interarrivaltime of frames
else
    lambda = 1/mean(diff(timeOfFrameDownload));
    tlambda = generateTraffic(arrival, duration, cv, lambda,ar);
end
cst = ceil(cumsum(tlambda));
bwPerSecond = histc(cst,unique(cst));
% bwPerSecond = tlambda;

% Sij = repmat(Sij,1,fps);
% Sij=reshape(Sij',1,numel(Sij))'*framesize;
% 
% packetnumbers=0:length(bwPerSecond);
% downloadedData=[0 cumsum(bwPerSecond)'];
% lastframe = find(sum(bwPerSecond)>cumsum(Sij),1,'last');
% replayedFrameData = cumsum(Sij(1:lastframe));
% 
% % for i=1:length(Sij)
% for i = 1:lastframe
% %     yi = i * framesize;
%     timeOfFrameDownload(i) = findX(packetnumbers,downloadedData,replayedFrameData(i));
% end

% pps = 10;
% pduration = 1/pps;
% bwPerSecond = repmat(bwPerSecond,1,pps);
% bwPerSecond=reshape(bwPerSecond',1,numel(bwPerSecond))'*pduration;
framesize2 = framesize;
dmax = d * mu / framesize; % packets
% bwPerSecond = diff(replayedFrameData);
% mu = mu * framesize;
TPolicyduration = dmax / lambda; % setting TPolicyduration to average stalling length for D-policy
% tlambda = ones(1,duration)*pduration;

% policies:
% 1 - n-policy
% 2 - D-policy
% 3 - T-policy
received = 1;
replaytime = zeros(1,duration);
% replayed = 1;
buffer = zeros(1,duration);
time = zeros(1,duration);
buffertime = zeros(1,duration);
autocorrelationvalue = 0.89;
% autocorrelationlag = 1;


playing = zeros(1,duration);
% if strcmp(arrival,'M')
%     tlambda = exprnd(1/lambda,1,duration); % s
% elseif strcmp(arrival,'G')
%     m=1/lambda;
%     v = cv * cv * m * m;
%     mx = log((m^2)/sqrt(v+m^2));
%     sigma = sqrt(log(v/(m^2)+1));
%     tlambda = lognrnd(mx,sigma,1,duration);
% elseif strcmp(arrival,'AR')
%     m=1/lambda;
%     v = cv * cv * m * m;
%     mx = log((m^2)/sqrt(v+m^2));
%     sigma = sqrt(log(v/(m^2)+1));
%     tlambda = lognrnd(mx,sigma,1,duration);
%     for i = 2:duration
%         tlambda(i) = tlambda(i-1) * autocorrelationvalue + tlambda(i) * (1-autocorrelationvalue);
%     end
% end

if plotme
    
    windowSizea = 100;
    a1 = (1/windowSizea)*ones(1,windowSizea);
    figure(2);clf
    subplot(4,1,1);
    plot(cumsum(tlambda),'.-');
    subplot(4,1,2);
    plot(cumsum(tlambda),filter(a1,1,1./tlambda),'.-r');
    ylabel('network traffic')
    subplot(4,1,3);
    windowSizeb = 24;
    b1 = (1/windowSizeb)*ones(1,windowSizeb);
    plot((1:duration)/fps,filter(b1,1,tmu),'.-b');
    ylabel('video bit rate')
end

% if strcmp(service,'M')
%     tmu = exprnd(mu,1,duration); % s
% elseif strcmp(service,'G')
%     m=mu;
%     v = cv * cv * m * m;
%     mx = log((m^2)/sqrt(v+m^2));
%     sigma = sqrt(log(v/(m^2)+1));
%     tmu = lognrnd(mx,sigma,1,duration);
% elseif strcmp(service,'AR')
%     m=mu;
%     v = cv * cv * m * m;
%     mx = log((m^2)/sqrt(v+m^2));
%     sigma = sqrt(log(v/(m^2)+1));
%     tmu = lognrnd(mx,sigma,1,duration);
%     for i = 2:duration
%         tmu(i) = tmu(i-1) * autocorrelationvalue + tmu(i) * (1-autocorrelationvalue);
%     end
% end

frame = 1;
% tmu = exprnd(mu,1,duration); % packets per frame
% tmu = lognrnd(mx,sigma,1,duration);
% tmu = zeros(1,duration)+1/mu; % s
% tmu2 = tmu;
% tlambda = exprnd(1/lambda,1,duration); % s
% tlambda2 = tlambda;
% tlambda = zeros(1,duration)+1/lambda; % s
% tmu = (exprnd(1/mu,1,duration) + exprnd(1/mu,1,duration))/2; % s
% tlambda = lognrnd(1/lambda,(1),1,duration); % s
startstallingtime=0;
TP = false;

i = 1;
while received < length(bwPerSecond) && frame < length(Sij) && received < length(tlambda)
%while time <= 700
    i = i + 1;
    if playing(i-1) == 1
        % busy state
        if time(i - 1) + framesize2 < time(i - 1) + tlambda(received) && buffer(i-1) >= tmu(frame)
            % replay 1 frame
            time(i) = time(i - 1) + framesize2;
            buffer(i) = buffer(i - 1) - tmu(frame);
            buffertime(i) = max(0,buffertime(i - 1) - framesize);
            tlambda(received) = tlambda(received) - framesize2;
            frame = frame + 1;
            framesize2 = framesize;
            playing(i) = 1;
            replaytime(i) = 1;
            TP = false;
        else
            % receive 1 unit
            time(i) = time(i - 1) + tlambda(received);
            buffer(i) = buffer(i - 1) + bwPerSecond(received);
            
            tmuCS = tmu(frame);
            k = 0;
            while tmuCS < buffer(i) && frame + k < length(Sij)
                k = k + 1;
                %                 if frame + k > duration
                %                     break;
                %                 end
                tmuCS = tmuCS + tmu(frame + k);
            end
            buffertime(i)=k*framesize;
            
            %             tmuCS = cumsum(tmu(frame:end)); % optimize runtime
            %             buffertime(i) = sum(tmuCS<=buffer(i)) * framesize;
            framesize2 = max(0, framesize2 - tlambda(received));
            received = received + 1;
            playing(i) = 1;
        end
    else
        % idle state
        time(i) = time(i - 1) + tlambda(received);
        buffer(i) = buffer(i - 1) + bwPerSecond(received);
        
        tmuCS = tmu(frame);
        k = 0;
        while tmuCS < buffer(i) && frame + k < length(Sij)
            k = k + 1;
            %                 if frame + k > duration
            %                     break;
            %                 end
            tmuCS = tmuCS + tmu(frame + k);
        end
        buffertime(i)=k*framesize;
        %
        %         tmuCS = cumsum(tmu(frame:end));
        %         buffertime(i) = sum(tmuCS<=buffer(i)) * framesize;
        framesize2 = max(0, framesize2 - tlambda(received));
        received = received + 1;
        
    end
    if policy == 2
        if buffertime(i) >= d
            % D-policy: buffer full
            playing(i) = 1;
        end
    elseif policy == 1
        if buffer(i) >= dmax;
            % n-policy: buffer full
            playing(i) = 1;
        end
    elseif policy == 3
        if startstallingtime + TPolicyduration <= time(i) && buffer(i) >= tmu(frame)
            % T-policy: wait T seconds and check if buffer not empty
            playing(i) = 1;
        end
    end
    if buffer(i) <= dmin || buffer(i) < tmu(frame)
        % buffer empty
        playing(i) = 0;
        if TP == false
            startstallingtime=time(i);
        end
        TP = true;
    end
    
end

emptyfieldindex = find(time == 0);
emptyfieldindex(1) = [];
playing(emptyfieldindex)=[];
buffer(emptyfieldindex)=[];
buffertime(emptyfieldindex)=[];
replaytime(emptyfieldindex)=[];
time(emptyfieldindex)=[];
% tlambda(emptyfieldindex)=[];

if plotme
    subplot(4,1,4)
    %     figure(1);clf;
%     plot(time, buffer,'.-');
    hold all;
    plot(time, buffertime);
    if policy == 1
%         plot([time(1) max(time)],[dmax dmax]);
    elseif policy == 2
%         plot([time(1) max(time)],[d d]*mu);
    end
    %     plot([time(1) time(end)],[dmin dmin]);
    xlabel('time [s]')
    ylabel('buffer status')
    %         plot(time(playing==1),buffer(playing==1),'*r')
%     plot(time(replaytime==1),buffer(replaytime==1),'ok')
    % plot(time,time==tmu,'o')
    legend({'seconds'})
    xlim([0 max(time)]);
    subplot(4,1,3);
    xlim([0 max(time)]);
    subplot(4,1,2);
    xlim([0 max(time)]);
    subplot(4,1,1);
    xlim([0 max(time)]);
end


playing = [0 playing];
startplaying = time(diff(playing)==1);
stopplaying = time(diff(playing)==-1);
n = length(stopplaying);
if length(stopplaying) < length(startplaying)
    stopplaying(end+1) = time(end);
end

% lengthofPlayingEvent = sum(stopplaying - startplaying)/length(startplaying);
rate = n / time(end);
ratio = 1 - sum(stopplaying - startplaying) / time(end);
N = n/sum(stopplaying - startplaying);
videolength=sum(stopplaying - startplaying);

startplaying = time(diff(playing)==1);
stopplaying = time(diff(playing)==-1);
if length(stopplaying) == length(startplaying)
    startplaying(end+1) = time(end);
end

L = mean(startplaying(2:end) - stopplaying);
% T = sum(startplaying(2:end) - stopplaying) / videolength;

if isnan(L)
    L = 0;
end

StallingRatio = sum(startplaying(2:end) - stopplaying)/(sum((startplaying(2:end) - stopplaying))+videolength);

QoE = exp(-(0.15 * L + 0.19) * N * 30);
% QoE = exp(-0.15 * T - 0.19 * N);
% QoE = exp(-0.15 * T) + exp(-0.19 * N) - 1;

% rate = rate * 60; % 1/s -> 1/min

[E1vod,E2vod] = usereng(rate*60,ratio*100,'vod');
[E1live,E2live] = usereng(rate*60,ratio*100,'live');

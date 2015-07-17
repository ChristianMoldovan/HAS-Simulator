function [QoE,E1vod,E2vod,E1live,E2live,N,L,rate,ratio,StallingRatio] = runSim2(lambda,mu,d,cv,plotme,arrival,service,policy)
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

timer = 600;
duration = 20000;
framesize = 1/24; % s
framesize2 = framesize;
dmax = d * mu; % packets
mu = mu * framesize;
TPolicyduration = dmax / lambda; % setting TPolicyduration to average stalling length for D-policy
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
if strcmp(arrival,'M')
    tlambda = exprnd(1/lambda,1,duration); % s
elseif strcmp(arrival,'G')
    m=1/lambda;
    v = cv * cv * m * m;
    mx = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
    tlambda = lognrnd(mx,sigma,1,duration);
elseif strcmp(arrival,'AR')
    m=1/lambda;
    v = cv * cv * m * m;
    mx = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
    tlambda = lognrnd(mx,sigma,1,duration);
    for i = 2:duration
        tlambda(i) = tlambda(i-1) * autocorrelationvalue + tlambda(i) * (1-autocorrelationvalue);
    end
end

if plotme
    figure(2);clf
%     subplot(3,1,1);
%     plot(cumsum(tlambda),'.-');
%     subplot(4,1,2);
%     plot(tlambda,'.-r');
%     subplot(4,1,3);
%     plot(cumsum(tlambda),tlambda,'.-');
end

if strcmp(service,'M')
    tmu = exprnd(mu,1,duration); % s
elseif strcmp(service,'G')
    m=mu;
    v = cv * cv * m * m;
    mx = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
    tmu = lognrnd(mx,sigma,1,duration);
elseif strcmp(service,'AR')
    m=mu;
    v = cv * cv * m * m;
    mx = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
    tmu = lognrnd(mx,sigma,1,duration);
    for i = 2:duration
        tmu(i) = tmu(i-1) * autocorrelationvalue + tmu(i) * (1-autocorrelationvalue);
    end
end

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
% while i < duration
while i < duration
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
            buffer(i) = buffer(i - 1) + 1;
            
            tmuCS = tmu(frame);
            k = 0;
            while tmuCS < buffer(i) && frame + k < duration
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
        buffer(i) = buffer(i - 1) + 1;
        
        tmuCS = tmu(frame);
        k = 0;
        while tmuCS < buffer(i) && frame + k < duration
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

if plotme
    subplot(3,1,3)
    %     figure(1);clf;
%     plot(time, buffer,'.-');
    hold all;
    plot(time, buffertime);
    if policy == 1
        plot([time(1) time(end)],[dmax dmax]);
    elseif policy == 2
        plot([time(1) time(end)],[d d]);
    end
    %     plot([time(1) time(end)],[dmin dmin]);
    xlabel('time [s]')
    ylabel('buffer status [s]')
%             plot(time(playing==1),buffer(playing==1),'*r')
    plot(time(replaytime==1),buffertime(replaytime==1),'ok')
    % plot(time,time==tmu,'o')
%     legend({'"packets"','seconds'})

    xlim([0 time(end)]);
subplot(3,1,2)
    plot(time, buffer,'.-');
    hold all
            plot(time(playing==1),buffer(playing==1),'*r')
    xlim([0 time(end)]);
    xlabel('time [s]')
    ylabel('buffer status [kB]')
    subplot(3,1,1);
            plot(time,playing,'r')
    xlim([0 time(end)]);
    ylim([0 1.1])
%     saveas(gcf,'figs\sysmodel','eps2c');
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

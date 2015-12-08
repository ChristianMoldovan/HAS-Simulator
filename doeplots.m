clear
close all
results = {'MM_3','MG_3','GM_3','simtimeMM_3','simtimeMG_3','simtimeGM_3'};
% QoEagg = zeros(100,3,10,6);
% Nagg = zeros(100,3,10,6);
% Lagg = zeros(100,3,10,6);
% simtimeagg = zeros(100,3,10,6);
% realmuagg = zeros(100,3,10,6);
% reallambdaagg = zeros(100,3,10,6);
% realaagg = zeros(100,3,10,6);
% realcv1agg = zeros(100,3,10,6);
% realcv2agg = zeros(100,3,10,6);
% arrivalagg = zeros(100,3,10,6);
% serviceagg = zeros(100,3,10,6);
% 
% for i = 1:length(results)
%     load(['results/' results{i} '']);
%     QoEagg(:,:,:,i) = QoE;
%     Nagg(:,:,:,i) = N;
%     Lagg(:,:,:,i) = L;
%     simtimeagg(:,:,:,i) = simtime;
%     realmuagg(:,:,:,i) = realmu;
%     reallambdaagg(:,:,:,i) = reallambda;
%     realaagg(:,:,:,i) = reala;
%     realcv1agg(:,:,:,i) = realcv1;
%     realcv2agg(:,:,:,i) = realcv2;
% %     arrivalagg(:,:,:,i) = arrival;
% %     serviceagg(:,:,:,i) = service;
% end

QoEagg = zeros(1,100*3*10*3 + 100*3*8*3);
Nagg = zeros(1,100*3*10*3 + 100*3*8*3);
Lagg = zeros(1,100*3*10*3 + 100*3*8*3);
simtimeagg = zeros(1,100*3*10*3 + 100*3*8*3);
realmuagg = zeros(1,100*3*10*3 + 100*3*8*3);
reallambdaagg = zeros(1,100*3*10*3 + 100*3*8*3);
realaagg = zeros(1,100*3*10*3 + 100*3*8*3);
realcv1agg = zeros(1,100*3*10*3 + 100*3*8*3);
realcv2agg = zeros(1,100*3*10*3 + 100*3*8*3);

for i = 1:length(results)
    load(['results/' results{i} '']);
    insertpoint = 1 + (i-1)*100*3*10;
    insertlength = 100*3*10-1;
    if i > 3
        insertpoint = 1 + 100*3*10*3 + (i-4)*100*3*8;
        insertlength = 100*3*8-1;
    end
    QoEagg(insertpoint:insertpoint+insertlength) = QoE(:);
    Nagg(insertpoint:insertpoint+insertlength) = N(:);
    Lagg(insertpoint:insertpoint+insertlength) = L(:);
%     simtimeagg(insertpoint:insertpoint+insertlength) = simtime(:);
    realmuagg(insertpoint:insertpoint+insertlength) = realmu(:);
    reallambdaagg(insertpoint:insertpoint+insertlength) = reallambda(:);
    realaagg(insertpoint:insertpoint+insertlength) = reala(:);
    realcv1agg(insertpoint:insertpoint+insertlength) = realcv1(:);
    realcv2agg(insertpoint:insertpoint+insertlength) = realcv2(:);
    
    qoesize = size(QoE);
    repa = repmat(a,qoesize(1),1,qoesize(3));
    repv = repmat(v,qoesize(1),qoesize(2),1);
    if i > 3
        repv = repmat(v,qoesize(1),qoesize(2),qoesize(3));
    end
    repsimtime = repmat(simtime,qoesize(1),qoesize(2),1);
    if i < 4
        repsimtime = repmat(simtime,qoesize(1),qoesize(2),qoesize(3));
    end
    reparrival = repmat(arrival,qoesize(1),qoesize(2),qoesize(3));
    repservice = repmat(service,qoesize(1),qoesize(2),qoesize(3));
    aagg(insertpoint:insertpoint+insertlength) = repa(:);
    vagg(insertpoint:insertpoint+insertlength) = repv(:);
    simtimeagg(insertpoint:insertpoint+insertlength) = repsimtime(:);
    arrivalagg(insertpoint:insertpoint+insertlength) = reparrival(:);
    serviceagg(insertpoint:insertpoint+insertlength) = repservice(:);
    
%     simtimeagg(insertpoint:insertpoint+insertlength) = simtime(:);
%     arrivalagg{i} = arrival;
%     serviceagg{i} = service;
end



%% % maineffectsplot(QoEagg',{Nagg',Lagg',realmuagg',reallambdaagg',realcv1agg',realcv2agg'});
figure(1)
[a,b]=histc(simtimeagg,[0 40 150 600 inf]);
[c,d]=histc(realcv1agg,[0.5 2.1 5 inf]);
[e,f]=histc(realcv2agg,[0.5 2.1 5 inf]);
[g,h]=histc(realaagg,[0.8 0.9 1 1.1 1.2 inf]);


h=maineffectsplot(QoEagg',{h',b',vagg',d',f',serviceagg'},'varnames',{'offered load', 'simtime','input c_v','c_v1','c_v2','service process'});
% h=maineffectsplot(QoEagg',{aagg',b',vagg',d',f',serviceagg'},'varnames',{'offered load', 'simtime','input c_v','c_v1','c_v2','service process'});
% hold all
% ylim([1 5]);
% ylabel('mean MOS')
% figure(2)
% h=maineffectsplot(Lagg',{aagg',simtimeagg',vagg',arrivalagg',serviceagg'},'varnames',{'offered load', 'simtime','c_v','arrival process','service process'});
% ylabel('mean length of stalling events')
% figure(3)
% h=maineffectsplot(Nagg',{aagg',simtimeagg',vagg',arrivalagg',serviceagg'},'varnames',{'offered load', 'simtime','c_v','arrival process','service process'});
% ylabel('mean rate of stalling events')
% figure(2)
% figure(4)
% maineffectsplot(QoEagg',{aagg',simtimeagg',vagg',arrivalagg',serviceagg'},'varnames',{'offered load', 'simtime','c_v','arrival process','service process'},'statistics','std');
[yy,y]=find(aagg==0.8);
% maineffectsplot(QoEagg(y)',{b(y)',vagg(y)',arrivalagg(y)',serviceagg(y)'},'varnames',{'simtime','c_v','arrival process','service process'});
ylim([1 5]);
%%
uv = unique(vagg);
% for i=1:3
figure(1);clf;
    for i=1:length(uv)
        cdfplot(realcv2agg((vagg==uv(i)) & (serviceagg=='G') & (simtimeagg==600) & (aagg==1.2)));
        hold all
    end
%     cdfplot(realcv1agg(vagg==max(vagg)))
% end
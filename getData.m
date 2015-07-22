function [bwPerSecond,Sij]=getData(runID,bw_faktor)%%
n=350;
if nargin<2,bw_faktor=1;end
if nargin<1,runID=0;end

httpHeaderSizePerSegment=242; % byte

LayerSize_ij = textread('data/ToS_s2_seg_size.txt','',n); % size of SVC layer persegment (Byte); each segment corresponds to 2sec of video playtime
Sij = cumsum(LayerSize_ij,2)+httpHeaderSizePerSegment; % size of segment (Byte); each segment corresponds to 2sec of video playtime

tau=2; % duration of segment (sec)

m=700; % read only first 700 time instants and adjust
fileTrafficPattern = sprintf('./data/eval_10/statistics_proposed/%04d/recorded_limits.txt',runID);
bwPerSecond = textread(fileTrafficPattern,'',m); % kilobit/s, since start of DASH client

targetBW=sum(Sij(:,3)*8/1000)/sum(bwPerSecond);
bwPerSecond=bwPerSecond*targetBW; % we have now similiar bandwidth for all!

V = [0; cumsum(repmat(bwPerSecond,5,1)*1000/8)]*bw_faktor; % (byte)
timeV = (0:length(V)-1)';

%% Define Volume function
Vt = @(xi) interp1(timeV,V,xi); % total amount data V(t) received by client during time [0,t]
%pos=ginput(1); plot(pos(1),Vt(pos(1)),'rX')
Vtrev = @(yi) interp1(V,timeV,yi); %reverse function of V(t): Vtrev(volume)=time
%pos=ginput(1); plot(Vtrev(pos(2)),pos(2),'g*')

%%


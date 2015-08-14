
function timeOfFrameDownload = getTimeOfFrameDownload(bwPerSecond, Sij)

packetnumbers=0:length(bwPerSecond);
downloadedData=[0 cumsum(bwPerSecond)'];
lastframe = find(sum(bwPerSecond)>cumsum(Sij),1,'last');
replayedFrameData = cumsum(Sij(1:lastframe));

% for i=1:length(Sij)
for i = 1:lastframe
    %     yi = i * framesize;
    timeOfFrameDownload(i) = findX(packetnumbers,downloadedData,replayedFrameData(i));
end

end

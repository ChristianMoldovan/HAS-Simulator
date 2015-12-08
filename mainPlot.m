function h=mainPlot(H,group,str,ylab,varargin)

a=3.4;b=4.2;

alpha=0.05;
for i=1:length(group)
    h(i)=subplot(1,length(group),(i));

    edges=unique(group{i});
    m=zeros(1,length(edges));
    c=zeros(1,length(edges));
    for j=1:length(edges)
        ind=find(group{i}==edges(j));
        m(j)=nanmean(H(ind));
        c(j)=tinv(1-alpha/2,length(ind)-1)*nanstd(H(ind))/sqrt(length(ind));
    end
    a=min([m-c,a]);b=max([m+c,b]);
    errorbar(edges,m,c,'.-',varargin{:});
    if i==1,ylabel(ylab);end
    xlabel(str{i},'interpreter','none');
    hold all
end


for i=1:length(group)
    subplot(1,length(group),i);
    ylim([a b]);
%    set(gca,'ytick',a:0.2:b);
end

%% Display the objective and residual trajectories of one dataset.
function lossrecord=valueplot(datas,timeend,index)
lossdata=datas{1};
trdata=datas{2}(:);
residual=datas{3}(:);
data=[lossdata,residual];
if trdata(end)>=timeend
    i=find(trdata>=timeend,1,'first');
    if trdata(i)==timeend
        trdata=trdata(1:i);
        data=data(1:i,:);
    elseif i>1
        endpoint=interp1(trdata(i-1:i),data(i-1:i,:),timeend,'linear');
        trdata=[trdata(1:i-1);timeend];
        data=[data(1:i-1,:);endpoint];
    end
end
lossrecord=data(end,:);
plotplt(data(:,1:2),trdata,data(:,3),index);
end

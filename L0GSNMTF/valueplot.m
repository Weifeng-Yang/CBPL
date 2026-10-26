function valueplot(datas,timeend,index)
tr=datas{2}(:);loss=datas{1};residual=datas{3}(:);
if tr(end)>=timeend
    i=find(tr>=timeend,1,'first');
    if tr(i)==timeend
        tr=tr(1:i);loss=loss(1:i,:);residual=residual(1:i);
    elseif i>1
        endpoint=interp1(tr(i-1:i),[loss(i-1:i,:),residual(i-1:i)],timeend,'linear');
        tr=[tr(1:i-1);timeend];loss=[loss(1:i-1,:);endpoint(1:2)];residual=[residual(1:i-1);endpoint(3)];
    end
end
plotplt(loss,tr,residual,index);
end

%% Display the recorded trajectory up to the observation time.
function plt0=valueplot(datas,timeend,name)
tr=datas{2}(:);
loss=datas{1};
R=datas{3}(:);
[tr,ids]=unique(tr,'stable');
loss=loss(ids,:);
R=R(ids);
x=[tr(tr<timeend);timeend];
if length(tr)==1
    x=tr;
    loss=loss(1,:);
    R=R(1);
else
    loss=interp1(tr,loss,x,'linear');
    R=interp1(tr,R,x,'linear');
end
plt0=plotplt(loss,x,R,name);
end

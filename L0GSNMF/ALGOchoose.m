%% data={objective values, running times, stationarity residuals}.
function [data,varss]=ALGOchoose(var,ngmar,aa,maxiteropt,Lapk,alpha,epsilon,delta,timeend)
[varss,loss,tr,residual]=CBPL(var,ngmar,aa,maxiteropt,Lapk,alpha,epsilon,delta,timeend);
data={loss,tr,residual};
end

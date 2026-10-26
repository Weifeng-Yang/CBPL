%% Run CBPL and return its objective values, times and certificates.
function [data,vars]=ALGOchoose(var,ngmar,aa,maxiteropt,Lapk,alpha,mu,epsilon,delta,timeend)
[vars,loss,tr,R]=CBPL(var,ngmar,aa,maxiteropt,Lapk,alpha,mu,epsilon,delta,timeend);
data={loss,tr,R};
end

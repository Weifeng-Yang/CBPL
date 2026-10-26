%% Solve the graph-regularized sparse matrix tri-factorization problem.
function [data,vars]=ALGOchoose(core,var,ngmar,coreaa,aa,maxiteropt,Lapk,alpha,mu,epsilon,delta,timeend)
[core,var,loss,tr,residual]=CBPL(core,var,ngmar,coreaa,aa,maxiteropt,Lapk,alpha,mu,epsilon,delta,timeend);
data{1}=loss;
data{2}=tr;
data{3}=residual;
vars{1}=core;
vars{2}=var;
end


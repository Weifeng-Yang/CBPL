%% Controlled block proximal linearized method for graph-regularized NMF.
function [var,loss,timerun,residual]=CBPL(var,ngmar,aa,maxiteropt,Lapk,alpha,epsilon,delta,timeend)
active=find(sum(var{end}.^2,1)>0);
if aa(end)==0
    J=sum(ngmar(:).^2)/2; loss=[J,J]; timerun=0; residual=0; return;
end
num=length(var);
Lmax=norm(Lapk,1);
for j=1:num
    varze{j}=zeros(size(var{j}));
end
loss(1,:)=compute(var,num,ngmar,Lapk,alpha,epsilon,active);
timerun=0;
residual(1,1)=certificate(var,ngmar,aa,num,Lapk,alpha,epsilon,active);

t1=tic;
for i=1:maxiteropt
    for j=1:num
        [U,ell]=grad(var,ngmar,j,num,Lapk,alpha,epsilon,active);
        ck=safeCBPL(var,j,U,ell,num,Lmax,alpha,epsilon,delta,Lapk,aa(end));
        a=aa(j);
        if isinf(a), a=[]; end
        var{j}=PROX(varze{j},var{j}-U/ck,a);
    end
    loss(i+1,:)=compute(var,num,ngmar,Lapk,alpha,epsilon,active);
    residual(i+1,1)=certificate(var,ngmar,aa,num,Lapk,alpha,epsilon,active);
    timerun(i+1)=toc(t1);
    if loss(i+1,2)>loss(i,2)+1e-9*max(1,loss(i,2))
        error('CBPL:Descent','The controlled energy increased.');
    end
    if stopcheck(timerun,timeend)
        break;
    end
end
end

function R=certificate(var,ngmar,aa,num,Lapk,alpha,epsilon,active)
R=0;
for j=1:num
    U=grad(var,ngmar,j,num,Lapk,alpha,epsilon,active);
    a=aa(j);
    if isinf(a), a=[]; end
    r=normalCBPL(var{j},U,a);
    R=R+sum(r(:).^2);
end
R=sqrt(R);
end

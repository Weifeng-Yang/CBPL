%% Controlled block proximal linearized method for sparse nonnegative CP.
function [var,loss,timerun,R]=CBPL(var,ngmar,aa,maxiteropt,Lapk,alpha,mu,epsilon,delta,timeend)
active=find(sum(var{end}.^2,1)>0);
if aa(end)==0
    J=sum(ngmar(:).^2)/2; loss=[J,J]; timerun=0; R=0; return;
end
num=length(var);
Lmax=norm(Lapk,1);
loss=compute(var,ngmar,Lapk,alpha,mu,epsilon,active);
timerun=0;
R=certificate(var,ngmar,aa,Lapk,alpha,mu,epsilon,active);
t1=tic;

for iter=1:maxiteropt
    for j=1:num
        [U,ell]=grad(var,ngmar,j,Lapk,alpha,mu,epsilon,active);
        ck=safeCBPL(var,j,U,ell,Lmax,alpha,mu,epsilon,delta,Lapk,aa(end));
        a=aa(j);
        if isinf(a), a=[]; end
        var{j}=PROX(zeros(size(var{j})),var{j}-U/ck,a);
    end
    loss(iter+1,:)=compute(var,ngmar,Lapk,alpha,mu,epsilon,active);
    R(iter+1,1)=certificate(var,ngmar,aa,Lapk,alpha,mu,epsilon,active);
    timerun(iter+1)=toc(t1);
    if loss(iter+1,2)>loss(iter,2)+1e-9*max(1,loss(iter,2))
        error('CBPL:Descent','Controlled objective increased.');
    end
    if mod(iter,100)==0
        fprintf('CBPL k=%d t=%.3f J=%.6g J+C=%.6g R=%.3g\n',iter,timerun(end),loss(end,1),loss(end,2),R(end));
    end
    if stopcheck(timerun,timeend), break; end
end
end

function R=certificate(var,ngmar,aa,Lapk,alpha,mu,epsilon,active)
R=0;
for j=1:length(var)
    U=grad(var,ngmar,j,Lapk,alpha,mu,epsilon,active);
    a=aa(j);
    if isinf(a), a=[]; end
    r=normalCBPL(var{j},U,a);
    R=R+sum(r(:).^2);
end
R=sqrt(R);
end

%% Controlled block proximal linearization with one safe update per block.
function [core,var,loss,tr,residual]=CBPL(core,var,ngmar,coreaa,aa,maxiteropt,Lapk,alpha,mu,epsilon,delta,timeend)
active=find(sum(var{end}.^2,1)>0);
if aa(end)==0
    J=sum(ngmar(:).^2)/2; loss=[J,J]; tr=0; residual=0; return;
end
Lmax=norm(Lapk,1);
loss(1,:)=compute(core,var,ngmar,Lapk,alpha,mu,epsilon,active);
tr=0;
residual(1)=certificate(core,var,ngmar,coreaa,aa,Lapk,alpha,mu,epsilon,active);
timerun=tic;
for iter=1:maxiteropt
    for j=0:2
        [U,ell]=grad(core,var,ngmar,j,Lapk,alpha,mu,epsilon,active);
        ck=safeCBPL(core,var,j,U,ell,Lmax,alpha,mu,epsilon,delta,Lapk,aa(end));
        if j==0, A=core; a=coreaa; else, A=var{j}; a=aa(j); end
        if isinf(a), a=[]; end
        Z=PROX(zeros(size(A)),A-U/ck,a);
        if j==0, core=Z; else, var{j}=Z; end
    end
    loss(iter+1,:)=compute(core,var,ngmar,Lapk,alpha,mu,epsilon,active);
    residual(iter+1)=certificate(core,var,ngmar,coreaa,aa,Lapk,alpha,mu,epsilon,active);
    if loss(iter+1,2)>loss(iter,2)+1e-9*max(1,loss(iter,2))
        error('CBPL:Descent','Controlled objective increased.');
    end
    tr(iter+1)=toc(timerun);
    if mod(iter,100)==0
        fprintf('k=%d t=%.3f J=%.6g J+C=%.6g R=%.3g\n',iter,tr(end),loss(end,1),loss(end,2),residual(end));
    end
    if stopcheck(tr,timeend), break; end
end
end

function R=certificate(core,var,ngmar,coreaa,aa,Lapk,alpha,mu,epsilon,active)
R=0;
for j=0:2
    U=grad(core,var,ngmar,j,Lapk,alpha,mu,epsilon,active);
    if j==0, A=core; a=coreaa; else, A=var{j}; a=aa(j); end
    if isinf(a), a=[]; end
    r=normalCBPL(A,U,a);
    R=R+sum(r(:).^2);
end
R=sqrt(R);
end



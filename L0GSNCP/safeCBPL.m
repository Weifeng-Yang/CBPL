%% Explicit safe curvature for one CP block update.
function ck=safeCBPL(var,j,U,ell,Lmax,alpha,mu,epsilon,delta,Lapk,budget)
num=length(var);
if nargin<11, budget=numel(var{end}); end
if j==num
    A=var{end};
    if budget<numel(A)
        active=find(sum(A.^2,1)>0);
        if isempty(active), ck=2.02*delta; return; end
        p=project(A,Lapk); z=A-p;
        q=sum(p(:,active).^2,1); r=sqrt(q);
        zn=sqrt(sum(z(:,active).^2,1));
        Theta=[1/64,1/32,1/16,1/8,1/4,1/2];
        cs=inf;
        for theta=Theta
            a=(1-theta)^2*q; b=(zn+theta*r).^2;
            aa=6./a.^2+20*b./a.^3; dd=2./a.^2; bb=8*sqrt(b)./a.^2.5;
            hh=(aa+dd+hypot(aa-dd,2*bb))/2;
            t=theta*min(r);
            c=max([2*delta,2*norm(U,'fro')/t,ell+alpha*Lmax+delta+epsilon*max([hh+18/5,1])]);
            cs=min(cs,c);
        end
        ck=1.01*cs;
    else
    A=var{num};
    p=project(A,Lapk);
    z=A-p;
    q=sum(p.^2,1);
    r=sqrt(q);
    zn=sqrt(sum(z.^2,1));
    gn=sqrt(sum(U.^2,1));
    Theta=[1/64,1/32,1/16,1/8,1/4,1/2];
    cs=inf(size(q));
    for theta=Theta
        a=(1-theta)^2*q;
        b=(zn+theta*r).^2;
        aa=6./a.^2+20*b./a.^3;
        dd=2./a.^2;
        bb=8*sqrt(b)./a.^2.5;
        hh=(aa+dd+hypot(aa-dd,2*bb))/2;
        c=max([2*delta*ones(size(q));gn./(theta*r);ell+alpha*Lmax+delta+epsilon*(hh+18/5)],[],1);
        cs=min(cs,c);
    end
    ck=1.01*max(cs);
    end
else
    A=var{j};
    rho=max(1,norm(A,'fro'));
    ts=[1,rho/64,rho/16,rho/4,rho,4*rho];
    c=inf;
    for t=ts
        if num==2
            ctrl=0;
        elseif j<num-1
            s=sum(var{j}.^2,1);
            ref=sum(var{num-1}.^2,1);
            ctrl=mu*max(3*(sqrt(s)+t).^2-ref);
        else
            s=sum(var{j}.^2,1);
            ref=zeros(size(s));
            for k=1:num-2
                ref=ref+sum(var{k}.^2,1);
            end
            ctrl=mu*max(3*(num-2)*(sqrt(s)+t).^2-ref);
        end
        c=min(c,max([2*delta,2*norm(U,'fro')/t,ell+ctrl+delta]));
    end
    ck=1.01*c;
end
if ~isfinite(ck)||ck<=0, error('CBPL:Curvature','Invalid safe curvature.'); end
end

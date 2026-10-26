%% Block gradient and fitting curvature; block zero is the independent core.
function [U,ell]=grad(core,var,ngmar,j,Lapk,alpha,mu,epsilon,active)
if nargin<9, active=find(sum(var{end}.^2,1)>0); end
inactive=setdiff(1:size(var{end},2),active);
A=var{1}; B=var{2};
if j==0
    U=(A'*A)*core*(B'*B)-A'*(ngmar*B);
    ell=largest(A'*A)*largest(B'*B);
    d=sum(A.^2,1)-sum(core.^2,2)';
    U=U-mu*core.*d';
elseif j==1
    T=core*(B'*B)*core';
    U=A*T-(ngmar*B)*core'; ell=largest(T);
    d=sum(A.^2,1)-sum(core.^2,2)';
    U=U+mu*A.*d;
else
    T=core'*(A'*A)*core;
    U=B*T-(ngmar'*A)*core; ell=largest(T);
    p=project(B,Lapk);
    q=sum(p.^2,1); s=sum(B.^2,1); 
    if any(q(active)<=0), error('CBPL:Domain','Active graph columns must stay nonzero.'); end
    h=s(active)./q(active).^2;
    U=U+alpha*(Lapk*B);
    U(:,active)=U(:,active)+epsilon*(2*B(:,active)./q(active).^2-4*s(active).*p(:,active)./q(active).^3).*(1-1./h.^2);
end
if j==0, U(:,inactive)=U(:,inactive)+epsilon*core(:,inactive);
elseif j==2, U(:,inactive)=U(:,inactive)+epsilon*var{2}(:,inactive); end
end

function ell=largest(T)
ell=max(0,max(eig((T+T')/2)));
end

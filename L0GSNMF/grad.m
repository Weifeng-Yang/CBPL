%% Mixed-state block gradient and fitting curvature.
function [U,ell]=grad(var,ngmar,j,num,Lapk,alpha,epsilon,active)
if nargin<8, active=find(sum(var{end}.^2,1)>0); end
inactive=setdiff(1:size(var{end},2),active);
k=3-j;
T=var{k}'*var{k};
ell=max(0,max(eig((T+T')/2)));
if j==1
    U=var{1}*T-ngmar*var{2};
else
    U=var{2}*T-ngmar'*var{1};
end
if j==num
    A=var{num};
    p=project(A,Lapk);
    q=sum(p.^2,1);
    s=sum(A.^2,1);
    
    h=s(active)./q(active).^2;
    if any(q(active)<=0)
        error('CBPL:Domain','q must stay positive.');
    end
    U=U+alpha*(Lapk*A);
    U(:,active)=U(:,active)+epsilon*(2*A(:,active)./q(active).^2-4*s(active).*p(:,active)./q(active).^3).*(1-1./h.^2);
end
U(:,inactive)=U(:,inactive)+epsilon*var{j}(:,inactive);
end

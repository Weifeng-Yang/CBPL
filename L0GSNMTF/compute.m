%% Calculate the original objective and the controlled function.
function loss=compute(core,var,ngmar,Lapk,alpha,mu,epsilon,active)
if nargin<8, active=find(sum(var{end}.^2,1)>0); end
inactive=setdiff(1:size(var{end},2),active);
A=var{1}; B=var{2}; D=A'*(ngmar*B);
fit=(sum(ngmar(:).^2)+sum(sum((A'*A)*core*(B'*B).*core))-2*sum(D(:).*core(:)))/2;
d=sum(A.^2,1)-sum(core.^2,2)';
balance=mu/4*sum(d.^2);
LA=Lapk*B;
graph=alpha/2*sum(sum(B.*LA));
p=project(B,Lapk);
q=sum(p.^2,1); s=sum(B.^2,1);

if any(q(active)<=0), error('CBPL:Domain','Active graph columns must stay nonzero.'); end
h=s(active)./q(active).^2;
C=balance+epsilon*sum(h+1./h-2);
tail=sum(sum(var{end}(:,inactive).^2))+sum(sum(core(:,inactive).^2));
C=C+epsilon/2*tail;
J=fit+graph;
loss=[J,J+C];
end

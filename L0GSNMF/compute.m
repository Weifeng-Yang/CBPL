%% Original objective J and controlled energy J+C.
function loss=compute(var,num,ngmar,Lapk,alpha,epsilon,active)
if nargin<7, active=find(sum(var{end}.^2,1)>0); end
inactive=setdiff(1:size(var{end},2),active);
T=var{1}'*var{1};
D=var{2}'*var{2};
Z=ngmar*var{2};
fit=(sum(ngmar(:).^2)+sum(T(:).*D(:))-2*sum(sum(var{1}.*Z)))/2;
A=var{num};
LA=Lapk*A;
graph=alpha/2*sum(sum(A.*LA));
p=project(A,Lapk);
q=sum(p.^2,1);
s=sum(A.^2,1);

h=s(active)./q(active).^2;
if any(q(active)<=0)
    error('CBPL:Domain','The last projected columns must be nonzero.');
end
C=epsilon*sum(h+1./h-2);
tail=0;
for k=1:length(var), tail=tail+sum(sum(var{k}(:,inactive).^2)); end
C=C+epsilon/2*tail;
J=fit+graph;
loss=[J,J+C];
if any(~isfinite(loss))
    error('CBPL:Numeric','The objective value is nonfinite.');
end
end

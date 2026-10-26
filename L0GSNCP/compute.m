%% Calculate the original and controlled objective values.
function loss=compute(var,ngmar,Lapk,alpha,mu,epsilon,active)
if nargin<7, active=find(sum(var{end}.^2,1)>0); end
inactive=setdiff(1:size(var{end},2),active);
num=length(var);
T=ones(size(var{1},2));
V=var{num};
for k=num-1:-1:2
    A=var{k};
    V=reshape(reshape(A,size(A,1),1,[]).*reshape(V,1,size(V,1),[]),[],size(A,2));
end
for k=1:num
    T=T.*(var{k}'*var{k});
end
Z=reshape(ngmar,size(var{1},1),[])*V;
fit=(sum(ngmar(:).^2)+sum(T(:))-2*sum(sum(var{1}.*Z)))/2;
balance=0;
for k=1:num-2
    d=sum(var{k}.^2,1)-sum(var{num-1}.^2,1);
    balance=balance+mu/4*sum(d.^2);
end
A=var{num};
LA=Lapk*A;
graphvalue=alpha/2*sum(sum(A.*LA));
p=project(A,Lapk);
q=sum(p.^2,1);
s=sum(A.^2,1);

h=s(active)./q(active).^2;
if any(q(active)<=0), error('CBPL:Domain','Last projected columns must be nonzero.'); end
C=balance+epsilon*sum(h+1./h-2);
tail=0;
for k=1:length(var), tail=tail+sum(sum(var{k}(:,inactive).^2)); end
C=C+epsilon/2*tail;
J=fit+graphvalue;
loss=[J,J+C];
if any(~isfinite(loss)), error('CBPL:Numeric','Nonfinite objective value.'); end
end

%% One mixed-state CP block gradient and fitting curvature.
function [U,ell]=grad(var,ngmar,j,Lapk,alpha,mu,epsilon,active)
if nargin<8, active=find(sum(var{end}.^2,1)>0); end
inactive=setdiff(1:size(var{end},2),active);
num=length(var);
T=ones(size(var{1},2));
index=setdiff(1:num,j);
V=var{index(end)};
for k=length(index)-1:-1:1
    A=var{index(k)};
    V=reshape(reshape(A,size(A,1),1,[]).*reshape(V,1,size(V,1),[]),[],size(A,2));
end
for k=index
    T=T.*(var{k}'*var{k});
end
X=reshape(permute(ngmar,[j,index]),size(var{j},1),[]);
U=var{j}*T-X*V;
ell=max(0,max(eig((T+T')/2)));
if j<num
    if j<num-1
        d=sum(var{j}.^2,1)-sum(var{num-1}.^2,1);
    else
        d=zeros(1,size(var{j},2));
        for k=1:num-2
            d=d+sum(var{j}.^2,1)-sum(var{k}.^2,1);
        end
    end
    U=U+mu*var{j}.*d;
end
if j==num
    A=var{num};
    p=project(A,Lapk);
    q=sum(p.^2,1);
    s=sum(A.^2,1);
    
    h=s(active)./q(active).^2;
    if any(q(active)<=0), error('CBPL:Domain','q must stay positive.'); end
    U=U+alpha*(Lapk*A);
    U(:,active)=U(:,active)+epsilon*(2*A(:,active)./q(active).^2-4*s(active).*p(:,active)./q(active).^3).*(1-1./h.^2);
end
U(:,inactive)=U(:,inactive)+epsilon*var{j}(:,inactive);
end

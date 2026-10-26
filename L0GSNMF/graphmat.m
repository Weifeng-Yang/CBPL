%% Graph Laplacian
function Lapk=graphmat(ngmar)
dimension=size(ngmar); n=dimension(end); X=reshape(ngmar,[],n);
s=full(sum(X.^2,1)); k=min(5,n-1); ids=zeros(k,n); dd=zeros(k,n);
for a=1:250:n
    id=a:min(n,a+249); D=max(0,s'+s(id)-2*full(X'*X(:,id)));
    for j=1:length(id), D(id(j),j)=inf; end
    [dd(:,id),ids(:,id)]=mink(D,k,1);
end
scale=median(dd(dd>0)); if isempty(scale)||scale<=0, scale=1; end
W=sparse(repelem(1:n,k),ids(:)',exp(-dd(:)'/scale),n,n);
W=max(W,W');
degree=full(sum(W,2));
Lapk=(spdiags(degree,0,n,n)-W)/max(degree);
end

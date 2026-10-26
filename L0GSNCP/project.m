%% Orthogonal range(L) projection for every connected graph component.
function p=project(A,Lapk)
persistent Lsave labels counts
if isempty(Lsave)||~isequal(Lsave,Lapk)
    W=-Lapk;
    W=W-spdiags(diag(W),0,size(W,1),size(W,1));
    labels=conncomp(graph(W))';
    counts=accumarray(labels,1);
    Lsave=Lapk;
end
if length(counts)==1
    p=A-mean(A,1);
else
    p=A;
    for j=1:size(A,2)
        averages=accumarray(labels,A(:,j))./counts;
        p(:,j)=A(:,j)-averages(labels);
    end
end
end

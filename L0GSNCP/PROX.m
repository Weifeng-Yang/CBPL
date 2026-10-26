%% Nonnegative projection with an elementwise hard sparsity budget.
function [X,IX]=PROX(X,U,a)
S=max(U,0);
if isempty(a)
    X=S;
    IX=find(X>0);
    return;
end
X(:)=0;
[~,IX]=maxk(S(:),a);
X(IX)=S(IX);
end

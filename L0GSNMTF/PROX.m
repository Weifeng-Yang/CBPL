%% Nonnegative hard-thresholding projection.
function [X,IX]=PROX(X,U,a)
S=max(U,0);
if isempty(a)
    X=S;
    IX=find(X>0);
    return;
end
X(:)=0;
IX=zeros(0,1);
if a>0
    [~,IX]=maxk(S(:),a);
    X(IX)=S(IX);
end
end

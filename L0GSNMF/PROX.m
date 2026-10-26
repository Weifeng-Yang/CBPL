%% Original nonnegative element-l0 projection.
function [X,IX]=PROX(X,U,a)
X=double(X);
U=double(U);
if ~isequal(size(X),size(U)) || ~isreal(U) || any(~isfinite(U(:)))
    error('PROX:Input','U must be finite real and match X.');
end
S=max(U,0);
if isempty(a)
    X=S;
    IX=find(X>0);
    return;
end
if ~isnumeric(a) || ~isscalar(a) || ~isreal(a) || ~isfinite(a) || a~=fix(a) || a<0 || a>numel(U)
    error('PROX:Budget','The element budget is invalid.');
end
X(:)=0;
IX=zeros(0,1);
if a==0, return; end
[~,IX]=maxk(S(:),a);
X(IX)=S(IX);
end

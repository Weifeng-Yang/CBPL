%% Minimum-norm residual in U plus the original regular normal cone.
function [r,n]=normalCBPL(X,U,a)
X=double(X);
U=double(U);
if ~isequal(size(X),size(U)) || ~isreal(X) || ~isreal(U) || ...
        any(~isfinite(X(:))) || any(~isfinite(U(:))) || any(X(:)<0)
    error('normalCBPL:Input','X must be feasible and U finite with the same shape.');
end
r=U;
IZ=(X==0);
r(IZ)=min(U(IZ),0);
if ~isempty(a)
    if ~isnumeric(a) || ~isscalar(a) || ~isreal(a) || ~isfinite(a) || a~=fix(a) || a<0 || a>numel(X)
        error('normalCBPL:Budget','The element budget is invalid.');
    end
    active=nnz(X);
    if active>a
        error('normalCBPL:Feasibility','X exceeds its element budget.');
    end
    if active==a
        r(IZ)=0;
    end
end
n=r-U;
if any(~isfinite(r(:))) || any(~isfinite(n(:)))
    error('normalCBPL:Numeric','The normal residual is nonfinite.');
end
end

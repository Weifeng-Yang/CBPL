%% Minimum-norm residual in U plus the regular normal of the original set.
function [r,n]=normalCBPL(X,U,a)
r=U;
IZ=(X==0);
r(IZ)=min(U(IZ),0);
if ~isempty(a) && nnz(X)==a
    r(IZ)=0;
end
n=r-U;
end

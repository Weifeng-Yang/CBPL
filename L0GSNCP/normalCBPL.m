% Minimum-norm residual in U plus the true Frechet normal of the original set.
function [r,n]=normalCBPL(X,U,a,groups)
% r is the minimum-norm element of U+Nhat_E(X); n=r-U belongs to Nhat_E(X).
% Slack patterns use the regular normal.
if nargin<4
    groups=[];
end
X=double(X);
U=double(U);
if ~isequal(size(X),size(U)) || ~isreal(X) || ~isreal(U) || ...
        any(~isfinite(X(:))) || any(~isfinite(U(:))) || any(X(:)<0)
    error('normalCBPL:Feasibility','X must be finite nonnegative and U finite with the same shape.');
end
r=U;
IZ=(X==0);
r(IZ)=min(U(IZ),0);
if isempty(a)
    n=r-U;
    return;
end
if isempty(groups)
    ngroups=numel(X);
else
    groups=groupsCBPL(groups,numel(X));
    ngroups=numel(groups);
end
if ~isnumeric(a) || ~isscalar(a) || ~isreal(a) || ~isfinite(a) || a~=fix(a) || a<0 || a>ngroups
    error('normalCBPL:Budget','Invalid sparsity budget.');
end
if isempty(groups)
    active=nnz(X);
    if active>a
        error('normalCBPL:Feasibility','X exceeds its element budget.');
    end
    if active==a
        r(IZ)=0;
    end
else
    active=false(ngroups,1);
    for i=1:ngroups
        active(i)=any(X(groups{i})>0);
    end
    if nnz(active)>a
        error('normalCBPL:Feasibility','X exceeds its group budget.');
    end
    if nnz(active)==a
        for i=1:ngroups
            if ~active(i)
                r(groups{i})=0;
            end
        end
    end
end
n=r-U;
if any(~isfinite(r(:))) || any(~isfinite(n(:)))
    error('normalCBPL:Overflow','The normal residual overflowed.');
end
end

function groups=groupsCBPL(groups,n)
if ~iscell(groups)
    error('normalCBPL:Groups','groups must be a cell array of disjoint coordinate sets.');
end
used=false(n,1);
for i=1:numel(groups)
    ids=groups{i};
    if ~isnumeric(ids) || ~isvector(ids) || isempty(ids) || ~isreal(ids) || ...
            any(~isfinite(ids(:))) || any(ids(:)~=fix(ids(:))) || any(ids(:)<1) || any(ids(:)>n)
        error('normalCBPL:Groups','Each group must contain finite valid MATLAB linear indices.');
    end
    ids=ids(:);
    if numel(unique(ids))~=numel(ids) || any(used(ids))
        error('normalCBPL:Groups','Groups must be disjoint and contain no repeated coordinates.');
    end
    used(ids)=true;
    groups{i}=ids;
end
if ~all(used)
    error('normalCBPL:Groups','Groups must partition all coordinates.');
end
end

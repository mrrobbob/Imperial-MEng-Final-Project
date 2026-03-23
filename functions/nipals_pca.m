function [T,P,expVar] = nipals_pca(X,LV)
% [T,P,expVar] = nipals_pca(X,LV)

% Carl Emil Eskildsen, 2018
% handles missing values

% References
% Wold, et al., Principal Component Analysis, Chemometrics and intelligent Laboratory Systems, 1987, 2, 37-52
% Smilde, A., R. Bro and P. GELADI. 2004. Two‐way component and regression models. Pages 35‐56 in Multi‐Way Analysis with Applications in the Chemical Sciences. 1st ed. A. Smilde, R. Bro and P. GELADI eds. John Wiley & Sons, Ltd, West Sussex, UK.
% Bro & Smilde, Principal Component Analysis, AnalyticalMethods, 6, 2014, 2812

M = isnan(X);
X(M) = 0;


SSt = sum(diag(X'*X));      % total sum of squares

[n,m] = size(X);
itmax = 10000;              % max iterations
tol = 1e-4;                 % tolerance

T = zeros(n,LV);            % preallocate matrix for scores
P = zeros(m,LV);            % preallocate matrix for loadings
expVar = zeros(LV,2);       % prealocate vector for explained variance

for i = 1:LV
    it = 0;                 % count iterations
    [~,j] = max(diag(X'*X));
    tit = X(:,j);
    sse = tol+1;
    
    while sse > tol^2
        it = it+1;
        T(:,i) = tit;
        P(:,i) = X'*tit;
        pp = sqrt(P(:,i)'*P(:,i));
        P(:,i) = P(:,i)/pp;
        tit = X*P(:,i);
        sse = (T(:,i)-tit)'*(T(:,i)-tit);
        if itmax < it
            T(:,i) = tit;
            Mdisp('convergence has not been reached')
            break
        end
    end
    
    X = X-T(:,i)*P(:,i)';
    while X(M)'*X(M) > 1e-10
        X(M) = 0;
        Pt = T(:,i)*inv(T(:,i)'*T(:,i))*T(:,i)';
        Pp = P(:,i)*inv(P(:,i)'*P(:,i))*P(:,i)';
        X = (eye(n)-Pt)*X*(eye(m)-Pp);
    end
    SSe = sum(diag(X'*X));
    expVar(i,1) = (1-SSe/SSt)*100;
    if i==1
        expVar(i,2) = expVar(i,1);
    else
        expVar(i,2) = expVar(i,1)-expVar(i-1,1);
    end
end

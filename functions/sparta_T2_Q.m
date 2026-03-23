function [T2,Q] = T2_Q(X,T,P,LV)
% function [T2,Q] = T2_Q(X,P)
% Calculate Q-resuduals and Hotelling's T-squared
% INPUT     X:  centered data matrix
%           T:  scores
%           P:  loadings
%           LV: number of latent variables
%
% OUTPUT    T2: Hotelling's T-squared
%           Q:  Q-residuals

% Carl Emil Eskildsen, 2018

% Reference
% Momikos & MacGregor, Multivariate SPC Charts for Monitoring Batch Processes, Technometrics, 1995, 37:1, 41-59

[~,m] = size(X);
Q = diag(X*(eye(m)-P(:,1:LV)*P(:,1:LV)')*X'); % Q-residuals
% Q = diag(X*(eye(m)-P*P')*X'); % Q-residuals

Tnew = X*P(:,1:LV);
% [~,S,~] = svd(X);
% l = power(diag(S),2);
% T2 = diag(T(:,1:LV).*pinv(l(1:LV))*T(:,1:LV)'); % Hotelling's T-squared
T2 = diag(Tnew*pinv(T(:,1:LV)'*T(:,1:LV))*Tnew');

end
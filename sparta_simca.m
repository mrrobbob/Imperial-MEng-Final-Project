clear; clc; close all;

mncn = @(x) (x-mean(x));


J = 200;
I = 600;
S = nan(J,3);
sigma = 8;
mu = [40;100;160];
for i = 1:3
    S(:,i) = normpdf(1:J,mu(i),sigma)';
end
% plot(1:J,S); shg

C = zeros(I,3);
idx = [1;101;201;301;401;501];
j = [1;2;3;1;2;3];
for i = 1:size(idx,1)
    C(idx(i):idx(i)+99,j(i))=rand(100,1);
end
clear idx j sigma mu i

X = C*S'+1E-3.*rand(I,J);
% plot(1:J,mncn(X(1:100,:)))
%%
Xcal = X(1:300,:);
Xtest = X(301:end,:);

Xcal_mncn = mncn(Xcal(101:200,:));
Xtest_mncn = Xtest(1:100,:)-mean(Xcal(101:200,:));

LV = 1;
[Tcal,P,expVar] = nipals_pca(Xcal_mncn,LV);
[~,Q_cal] = T2_Q(Xcal_mncn,Tcal,P,LV);
[~,Q_test] = T2_Q(Xtest_mncn,Tcal,P,LV);
CI_UpperLim = sparta_upperCI_Qstat(Xcal_mncn,P);

figure;
hs1 = subplot(2,1,1); hold on;
hh1 = histogram(Q_cal);
line([CI_UpperLim,CI_UpperLim],ylim)

hs2 = subplot(2,1,2);
hh2 = histogram(Q_test);
linkaxes([hs1,hs2],'x');
set(hh2,'BinWidth',get(hh1,'BinWidth'));
line([CI_UpperLim,CI_UpperLim],ylim)
shg
%%
Ecal = Xcal_mncn-Tcal*P';

figure;
hs1 = subplot(2,1,1); hold on;
hh1 = histogram(Q_cal);

hs2 = subplot(2,1,2);
hh2 = histogram(Q_test);
linkaxes([hs1,hs2],'x');
set(hh2,'BinWidth',get(hh1,'BinWidth'));
shg

%%
Ecal = Xcal_mncn-Tcal*P';
lampda = eig(cov(Ecal));

theta = nan(3,1);
for i = 1:3
    theta(i) = sum(power(lampda,i));
end

g = theta(2)/theta(1);
h = power(theta(1),2)/theta(2);
% alpha = 0.05;
% z = norminv(1-alpha);
z = 1.96;
h0 = 1-(2*theta(1)*theta(3))/(3*power(theta(2),2));

beta = nan(2,1);
beta(1,1) = 1-theta(2)*h0*(1-h0)/power(theta(1),2); 
beta(2,1) = z*power(2*theta(2)*power(h0,2),0.5)/theta(1);
Qa = theta(1)*power(beta(1)+beta(2),1/h0);

Qaa = g*h*power(1-2/(9*h)+power(z*(2/9*h),0.5),3);



%%

% E = Proj*Xbg_mn(:,idx);
% %E = Xbg_mn;
% % qstat_bg = diag(E'*E);
% eigen = sort(eig(cov(E')),'descend');
% eigen = eigen(1:min(I,J)-1);
% 
% theta = nan(3,1);
% for i = 1:3
%     theta(i) = sum(power(eigen,i));
% end
% h0 = 1-(2*theta(1)*theta(3))/(3*power(theta(2),2));
% alpha = 0.05;
% beta = nan(2,1);
% beta(1,1) = norminv(1-alpha)*sqrt(2*theta(2)*power(h0,2))/theta(1);
% beta(2,1) = theta(2)*h0*(h0-1)/power(theta(1),2);
% Qa = theta(1)*power(beta(1)+1+beta(2),1/h0);
% upper_lim = Qa;
% % phat = NaN(1,1);
% 
% 
% figure;
% histogram(qstat_bg);
% line([Qa, Qa],ylim);

% %%
% % E = X-T(:,1:K)*P(:,1:K)';
% % E = Xbg(:,idx)'*Proj;
% % E = Proj*Xbg(:,idx);
% eigen = sort(eig(cov(E)),'descend');
% eigen = eigen(1:min(I,J)-2);
% theta = nan(3,1);
% for i = 1:3
%     theta(i) = sum(power(eigen,i));
% end
% 
% h0 = 1-(2*theta(1)*theta(3))/(3*power(theta(2),2));
% %h0 = 1-2/3;
% alpha = 0.05;
% beta = nan(2,1);
% beta(1,1) = norminv(1-alpha)*sqrt(2*theta(2)*power(h0,2))/theta(1);
% beta(2,1) = theta(2)*h0*(h0-1)/power(theta(1),2);
% Qa = theta(1)*power(beta(1)+1+beta(2),1/h0);
% 
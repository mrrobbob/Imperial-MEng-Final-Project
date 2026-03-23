clear; close all; clc;

pat_fun = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\code\functions';
pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S19P4WaterBck\short_acquisitions';

% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-01\S07P1PBSBckgr\short_acquisitions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-01\S10P1PBSBckgr\short_acquisitions';

mncn = @(x) (x-mean(x));
auto = @(x) (x-mean(x))./std(x);

cd(pat_fun);
load spartaalphawn.mat
[Xbck,ID,I,it] = sparta_import_shortacq(pat_data);
% [Xbck,~]=AutoSpike_Matrix(Xbck,wn,1);
Xbckp = sparta_preprocess(Xbck,wn);
% plot(wn,Xbckp);
%%

idx = false(I*it,1);
idx(randperm(I*it,ceil(I*it*0.2))) = true;

mnX = mean(Xbckp(idx,:));
E = Xbckp(~idx,:)-mnX;

uCI = sparta_upperCI_Qstat(E);
Q = diag(E*E');

histogram(Q);
line([uCI,uCI],ylim);
shg
clear E I ID idx it Q 
%%
figure;
pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S19P6MmixStoc\short_acquisitions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-01\S07P3Ljon0500\short_acquisitions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-01\S10P3SyneStoc0250\short_acquisitions';
cd(pat_fun);
[X,ID,I,it] = sparta_import_shortacq(pat_data);
% [X,~]=AutoSpike_Matrix(X,wn,1);
Xp = sparta_preprocess(X,wn);

% plot(wn,XL500p);

E = Xp-mnX;
Q = diag(E*E');
histogram(Q);
line([uCI,uCI],ylim);
shg
%%
idx = false(I*it,1);
idx(Q>uCI) = true;

plot(wn,Xp(idx,:))

%%

[~,a1] = min(abs(wn-2050));
[~,b1] = min(abs(wn-2150));

LV = 10;
[T1,P1,expVar1] = nipals_pca(mncn(Xp(idx,a1:b1)),LV);

[~,a2] = min(abs(wn-2150));
[~,b2] = min(abs(wn-2300));

[T2,P2,expVar2] = nipals_pca(mncn(Xp(idx,a2:b2)),LV);

%%
subplot(2,1,1);
plot(wn(a1:b1),P1(:,1)); shg

subplot(2,1,2);
plot(wn(a2:b2),P2(:,1:2));

%%
figure;
plot(1:sum(idx),auto(T1(:,1)),'.','markersize',10); hold on;
plot(1:sum(idx),auto(T2(:,1)),'.','markersize',10);
plot(1:sum(idx),auto(-T2(:,2)),'.','MarkerSize',10);
shg;
%%
xtemp = X(idx,:);

%%
figure;
for i = 250:280
    plot(wn,xtemp(i,:)); 
    pause;
end
    

%%
figure;
subplot(3,1,1);
plot(wn,Xbckp);
title('PBS');

subplot(3,1,2);
plot(wn,Xp);
title('samples');

subplot(3,1,3);
plot(wn,mean(Xbckp)); hold on;
plot(wn,mean(Xp),'--');
legend('PBS','Samples')
title('Mean spectra');
clear; close all; clc
pat_fun = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\code\functions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\check';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S18P3Syne0500\short_acquisitions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S16P5Syne1000\short_acquisitions';

% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-05-24\New_Tags_Jono_Double\short_acquisitions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-01\S09P8SyneStoc\short_acquisitions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-01\S14P5MrepStoc_test\short_acquisitions';

pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S15P2Ljon1000\short_acquisitions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S15P3Mrep1000\short_acquisitions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S16P5Syne1000\short_acquisitions';

% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-07\S26P6Mmix1000\short_acquisitions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-08\S32P6Mmix1000\short_acquisitions';

mncn = @(x) (x-mean(x));

cd(pat_fun);
load spartaalphawn.mat
[X,ID, I, it] = sparta_import_shortacq(pat_data);
% [X,~]=AutoSpike_Matrix(X,wn,1);

lambda_x = 1e4;     % smoothing parameter (generally 1e5 to 1e8) 
p_x = 0.0001;        % asymmetry parameter (generally 0.001)
d_x = 2;            % order of differences in penalty (generally 2)
Xbck = zeros(size(X));
for i = 1:size(X,1)
    Xbck(i,:) = [asysm(X(i,:)', lambda_x, p_x, d_x)]'; % Asymetric least squares to estimate baseline
%     plot(wn,X(i,:)-Xbck(i,:)); hold on;
%     plot(wn,X(i,:));
%     plot(wn,Xbck(i,:)); hold off;
%     pause
end

pol = 2;
win = 11;
der = 0;
Xp = sparta_SavGol_filter(X-Xbck,pol,win,der);

% plot all data
figure;
subplot(3,1,1);
plot(wn, X');
not_so_tight('Y'); grid on;

subplot(3,1,2);
plot(wn,Xp');
not_so_tight('Y'); grid on;
shg

% make PCA
LV = 1;
[T,P,expvar] = nipals_pca(mncn(Xp),LV);

subplot(3,1,3);
% mycolor_scatter([1:size(T,1)']',T*-1,repmat([1:it]',I,1));
mycolor_scatter([1:size(T,1)']',T,repmat([1:it]',I,1));
ylabel(['Score, LV1 (',num2str(expvar(1,1),'%.2f'),'%)']);
hcb = colorbar;
caxis([0 it]);
ylabel(hcb, ['Time [sec]'])
not_so_tight, grid on;


%%
I1 = 301;

[~,idx1] = min(abs(wn-2202));

figure;
subplot(3,1,1);
plot(1:it,mncn(X(I1:I1+it-1,idx1)),'.'); hold on;
plot(1:it,mncn(Xp(I1:I1+it-1,idx1)),'.'); hold off;

hs2 = subplot(3,1,2); hold on;
hs3 = subplot(3,1,3); hold on;
for i=I1:I1+it-1
    plot(hs2,wn,X(i,:)');
    plot(hs3,wn,Xp(i,:));
%     set(hs3,'Ylim',[0 800]);
    pause;
end

%%
a1 = 4;
a2 = 11;
a3 = 15;
figure;
subplot(2,1,1);
plot(wn,X([I1+a1-1,I1+a2-1,I1+a3-1],:));
subplot(2,1,2);
plot(wn,Xp([I1+a1-1,I1+a2,I1+a3-1],:))


%%

figure;
for i = 1:size(X,1)
    plot(wn,Xp(i,:)); grid on;
    line([2200, 2200],ylim); % Jon
    line([2216,2216], ylim); % yne
    line([2100,2100], ylim); % rep
    pause;
end
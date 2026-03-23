% figure_outline
clear; close all; clc

path_fun = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\code\functions';
path_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-08\S28P2Lrep1000\short_acquisitions';
path_figures = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\figures_for_report';

mncn = @(x) x-mean(x);

cd(path_fun);
load spartaalphawn.mat
[X,~, I, it] = sparta_import_shortacq(path_data);

% for i =1:100
%     plot(wn,X(i,:));
%     title(num2str(i));
%     pause;
% end
%% Figure, preprocessing
xtemp_spike = X(10,:);
xtemp_spike(1670) = xtemp_spike(1670)+1E3;
[xtemp,~]=AutoSpike_Matrix(xtemp_spike,wn,1);


lambda_x = 1e4;     % smoothing parameter (generally 1e5 to 1e8) 
p_x = 0.0001;        % asymmetry parameter (generally 0.001)
d_x = 2;            % order of differences in penalty (generally 2)
xbck = [asysm(xtemp', lambda_x, p_x, d_x)]'; % Asymetric least squares to estimate baseline

pol = 2;
win = 11;
der = 0;
xp = sparta_SavGol_filter(xtemp-xbck,pol,win,der);


figure;
subplot(4,1,1);
plot(wn,xtemp_spike);
grid on;

subplot(4,1,2);
plot(wn, xtemp); hold on;
plot(wn,xbck); hold off

subplot(4,1,3);
plot(wn, xtemp-xbck);

subplot(4,1,4);
plot(wn, xp);
xlabel('Raman shift (cm^-^1)')


% mysetsize(gcf,10,15);
% print(gcf,'-dtiff','-r1000',[path_figures,'\preprocessing.tiff']);
% cd(path_fun);
clear d_x der lambda_x p_x pol win xbck xp xtemp xtemp_spike I it X 
%% PCA plot
cd(path_fun);
path_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-08\S28P3Lrep0500\short_acquisitions';
[X,~, I, it] = sparta_import_shortacq(path_data);
Xp = sparta_preprocess(X,wn);

% make PCA
LV = 1;
[T,P,expvar] = nipals_pca(mncn(Xp),LV);

figure;
hs1 = subplot(2,1,1);
plot(wn,Xp); grid on;
not_so_tight('Y');
ylabel('Intensity');

hs2 = subplot(2,1,2);
% mycolor_scatter([1:size(T,1)']',T*-1,repmat([1:it]',I,1));
mycolor_scatter([1:size(T,1)']',T,repmat([1:it]',I,1),'off');
ylabel(['Score, LV1 (',num2str(expvar(1,1),'%.2f'),'%)']);
hcb = colorbar;
set(hcb,'Location','southoutside');
caxis([0 it]);
ylabel(hcb, ['Time [sec]'])
not_so_tight, grid on; shg

pos = get(hs2,'Position');
pos(2) = pos(2)+0.05;
set(hs2,'Position',pos);
mysetsize(gcf,15,20);

%% Visualisation of traps

figure;
subplot(4,1,1);
mycolor_scatter([541:600']',T(541:600),repmat([1:it]',3,1));
xlabel('Time (sec)')

hs1 = subplot(4,1,2);
mycolor(wn,Xp(541:560,:),[1:20]');
ylim([0 1300]);

hs2 = subplot(4,1,3);
mycolor(wn,Xp(561:580,:),[1:20]');
ylim([0 1300]);

hs3 = subplot(4,1,4);
mycolor(wn,Xp(581:600,:),[1:20]');
ylim([0 1300]);
pos = get(hs3,'Positio');
pause(0.1);

hcb = colorbar;
set(hcb,'Location','southoutside');
ylabel(hcb, ['Time [sec]'])
caxis([0 it]);
pause(0.1);
set(hs3,'Position',pos);
mysetsize(gcf,15,20);
%% Chi square - particle detection

cd(path_fun);
path_bck = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-08\S28P1WaterBck\short_acquisitions';
[Xbck,~, ~, ~] = sparta_import_shortacq(path_bck); % corresponding data with pure solvent
Xbckp = sparta_preprocess(Xbck,wn);

%% chi squared distribution, find particles
idx_bck = false(size(Xbckp,1),1);
idx_bck(randperm(size(Xbckp,1),ceil(size(Xbckp,1)*0.2))) = true;
mnX = mean(Xbckp(idx_bck,:)); % mean of background
E = Xbckp(~idx_bck,:)-mnX;
uCI = sparta_upperCI_Qstat(E); % this is the upper limit
Qbck = diag(E*E');
figure;
hh = histogram(Qbck); set(hh,'BinWidth',25000);
clear idx_bck E I
hl = line([uCI,uCI],ylim); shg
hl.LineWidth = 2;
xlim([2E5 8E5]);

Ex = Xp-mnX;
Q = diag(Ex*Ex');
idx = false(size(Xp,1),1);
idx(Q>uCI) = true;
text(4.5E5,25, ['Limit: ', num2str(uCI,'%.2e')]);


figure;
subplot(2,2,1);
plot(wn,Xp(21,:)); hold on;
plot(wn,mnX); hold off;
not_so_tight('Y'); grid on;
ylim([-50 250]);

subplot(2,2,2);
plot(wn,Xp(21,:)-mnX);
not_so_tight('Y'); grid on;
ylim([-50 250]);
text(1000,200, ['SSE: ', num2str(Q(21),'%.2e')]);


subplot(2,2,3);
plot(wn,Xp(40,:)); hold on;
plot(wn,mnX); hold off;
not_so_tight('Y'); grid on;
ylim([-50 250]);

subplot(2,2,4);
plot(wn,Xp(40,:)-mnX);
not_so_tight('Y'); grid on;
ylim([-50 250]);
text(1000,200, ['SSE: ', num2str(Q(40),'%.2e')]);

clear Ex expvar hcb hh hl hs1 hs2 hs3 idx it LV mnX P path_bck path_data pos Q Qbck T uCI X Xbck Xbckp Xp
%% plot probabilities of successful trap
cd(path_fun);
load M_matrix_20230602.mat
N_Traps02 = sparta_count_traps(M);
[T_Traps02,std02] = sparta_time_traps(M);
xlab = cell(size(ID_M,2),3);
for i = 1:size(xlab,1)
    xlab{i,1} = [ID_M{i}(6:9),ID_M{i}(10:13)];
end
clear ID_M M
idx = true(13,1);
idx(12) = false;
xlab = xlab(idx,:);
N_Traps02 = N_Traps02(idx);
T_Traps02 = T_Traps02(idx);
std02 = std02(idx);
clear idx

load M_matrix_20230607.mat
N_Traps07 = sparta_count_traps(M);
[T_Traps07,std07] = sparta_time_traps(M);
for i = 1:size(M,2)
    xlab{i,2} = [ID_M{i}(6:9),ID_M{i}(10:13)];
end
clear ID_M M

load M_matrix_20230608.mat
N_Traps08 = sparta_count_traps(M);
[T_Traps08,std08] = sparta_time_traps(M);
for i = 1:size(M,2)
    xlab{i,3} = [ID_M{i}(6:9),ID_M{i}(10:13)];
end
clear i M ID_M

[~,I] = sort(xlab(:,1));
N_Traps02 = N_Traps02(I);
T_Traps02 = T_Traps02(I);
std02 = std02(I);
[~,I] = sort(xlab(:,2));
N_Traps07 = N_Traps07(I);
T_Traps07 = T_Traps07(I);
std07 = std07(I);
[~,I] = sort(xlab(:,3));
N_Traps08 = N_Traps08(I);
T_Traps08 = T_Traps08(I);
std08 = std08(I);

xlab = xlab(I,3);
N_Traps = [N_Traps02;N_Traps07;N_Traps08];
T_Traps = [T_Traps02,T_Traps07,T_Traps08];
std_Traps = [std02;std07;std08]';
clear N_Traps02 N_Traps07 N_Traps08 I T_Traps02 T_Traps07 T_Traps08
N_Traps = N_Traps';

phat = nan(12,3);
CI = nan(12,6);
Idx = [1,3,5];
for i = 1:3
    [phat(:,i),CI(:,Idx(i):Idx(i)+1)] = sparta_CI_trap_prob(repmat(50,12,1),N_Traps(:,i));
end
CI(CI>1) = 1;

figure;
subplot(1,3,1);
errorbar([0.5,2.5,4.5],phat(1:3,1),phat(1:3,1)-CI(1:3,2),'s','linewidth',1,'MarkerSize',8); hold on;
errorbar([1,3,5],phat(1:3,2),phat(1:3,2)-CI(1:3,4),'s','linewidth',1,'MarkerSize',8); 
errorbar([1.5,3.5,5.5],phat(1:3,3),phat(1:3,3)-CI(1:3,6),'s','linewidth',1,'MarkerSize',8); 
grid on; set(gca,'XTick',[1,3,5],'XTickLabel',xlab(1:3),'XLim',[0, 6]);
ylim([0 1.1]);
axis square;
ylabel('Probability');

subplot(1,3,2);
errorbar([0.5,2.5,4.5],phat(4:6,1),phat(4:6,1)-CI(4:6,2),'s','linewidth',1,'MarkerSize',8); hold on;
errorbar([1,3,5],phat(4:6,2),phat(4:6,2)-CI(4:6,4),'s','linewidth',1,'MarkerSize',8); 
errorbar([1.5,3.5,5.5],phat(4:6,3),phat(4:6,3)-CI(4:6,6),'s','linewidth',1,'MarkerSize',8); 
grid on; set(gca,'XTick',[1,3,5],'XTickLabel',xlab(4:6),'XLim',[0, 6]);
ylim([0 1.1]);
axis square;
pos = get(gca,'Position');
pause(0.1);
legend('Day 1','Day 2','Day 3','Location','northoutside');
pause(0.1);
set(gca,'Position',pos);

subplot(1,3,3);
errorbar([0.5,2.5,4.5],phat(10:12,1),phat(10:12,1)-CI(10:12,2),'s','linewidth',1,'MarkerSize',8); hold on;
errorbar([1,3,5],phat(10:12,2),phat(10:12,2)-CI(10:12,4),'s','linewidth',1,'MarkerSize',8); 
errorbar([1.5,3.5,5.5],phat(10:12,3),phat(10:12,3)-CI(10:12,6),'s','linewidth',1,'MarkerSize',8); 
grid on; set(gca,'XTick',[1,3,5],'XTickLabel',xlab(10:12),'XLim',[0, 6]);
ylim([0 1.1]);
axis square;

mysetsize(gcf,10,18);

% clear CI i Idx N_Traps pos

%% average time in trap

figure;
subplot(1,3,1);
errorbar([0.5,2.5,4.5],T_Traps(1:3,1),std_Traps(1:3,1),'s','linewidth',1,'MarkerSize',8); hold on;
errorbar([1,3,5],T_Traps(1:3,2),std_Traps(1:3,2),'s','linewidth',1,'MarkerSize',8); 
errorbar([1.5,3.5,5.5],T_Traps(1:3,3),std_Traps(1:3,3),'s','linewidth',1,'MarkerSize',8); 
grid on; set(gca,'XTick',[1,3,5],'XTickLabel',xlab(1:3),'XLim',[0, 6]);
ylim([-1 11]);
axis square;
ylabel('Time (sec)');

subplot(1,3,2);
errorbar([0.5,2.5,4.5],T_Traps(4:6,1),std_Traps(4:6,1),'s','linewidth',1,'MarkerSize',8); hold on;
errorbar([1,3,5],T_Traps(4:6,2),std_Traps(4:6,2),'s','linewidth',1,'MarkerSize',8); 
errorbar([1.5,3.5,5.5],T_Traps(4:6,3),std_Traps(4:6,3),'s','linewidth',1,'MarkerSize',8); 
grid on; set(gca,'XTick',[1,3,5],'XTickLabel',xlab(4:6),'XLim',[0, 6]);
ylim([-1 11]);
axis square;
pos = get(gca,'Position');
pause(0.1);
legend('Day 1','Day 2','Day 3','Location','northoutside');
pause(0.1);
set(gca,'Position',pos);

subplot(1,3,3);
errorbar([0.5,2.5,4.5],T_Traps(10:12,1),std_Traps(10:12,1),'s','linewidth',1,'MarkerSize',8); hold on;
errorbar([1,3,5],T_Traps(10:12,2),std_Traps(10:12,2),'s','linewidth',1,'MarkerSize',8); 
errorbar([1.5,3.5,5.5],T_Traps(10:12,3),std_Traps(10:12,3),'s','linewidth',1,'MarkerSize',8); 
grid on; set(gca,'XTick',[1,3,5],'XTickLabel',xlab(10:12),'XLim',[0, 6]);
ylim([-1 11]);
axis square;

mysetsize(gcf,10,18);

%% Mixtures - probability of trap
load Mmix_Traps.mat
load SMLmix_Traps.mat
N_Traps = [N_Traps_Mmix1;N_Traps_Mmix2;N_Traps_Mmix3]';

phat = nan(3,3);
CI = nan(3,6);
Idx = [1,3,5];
for i = 1:3
    [phat(:,i),CI(:,Idx(i):Idx(i)+1)] = sparta_CI_trap_prob(repmat(50,3,1),N_Traps(:,i));
end
CI(CI>1) = 1;

figure;
subplot(1,2,1);
errorbar([0.5,2.5,4.5],phat(1:3,1),phat(1:3,1)-CI(1:3,2),'s','linewidth',1,'MarkerSize',8); hold on;
errorbar([1,3,5],phat(1:3,2),phat(1:3,2)-CI(1:3,4),'s','linewidth',1,'MarkerSize',8); 
errorbar([1.5,3.5,5.5],phat(1:3,3),phat(1:3,3)-CI(1:3,6),'s','linewidth',1,'MarkerSize',8); 
grid on; set(gca,'XTick',[1,3,5],'XTickLabel',ID,'XLim',[0, 6]);
ylim([0 1.1]);
axis square;
ylabel('Probability');
pos = get(gca,'Position');
pause(0.1);
legend('Day 1','Day 2','Day 3','Location','northoutside');
pause(0.1);
set(gca,'Position',pos);

n_tot = repmat(50,3,1);
n_suc = N_Traps(:,3);
p0 = mean(n_suc)/50;

[chi2_stat,df,p_value] = sparta_chi_square_test(n_suc,n_tot,p0);

%%

N_Traps = [N_Traps_SMLmix1([2,1,3]);N_Traps_SMLmix2;N_Traps_SMLmix3]';

phat = nan(3,3);
CI = nan(3,6);
Idx = [1,3,5];
for i = 1:3
    [phat(:,i),CI(:,Idx(i):Idx(i)+1)] = sparta_CI_trap_prob(repmat(50,3,1),N_Traps(:,i));
end
CI(CI>1) = 1;

subplot(1,2,2);
errorbar([0.5,2.5,4.5],phat(1:3,1),phat(1:3,1)-CI(1:3,2),'s','linewidth',1,'MarkerSize',8); hold on;
errorbar([1,3,5],phat(1:3,2),phat(1:3,2)-CI(1:3,4),'s','linewidth',1,'MarkerSize',8); 
errorbar([1.5,3.5,5.5],phat(1:3,3),phat(1:3,3)-CI(1:3,6),'s','linewidth',1,'MarkerSize',8); 
grid on; set(gca,'XTick',[1,3,5],'XTickLabel',{'L';'M';'S'},'XLim',[0, 6]);
ylim([0 1.1]);
axis square;
pos = get(gca,'Position');
pause(0.1);
pause(0.1);
set(gca,'Position',pos);

mysetsize(gcf,10,18);

n_tot = repmat(50,3,1);
n_suc = N_Traps(:,3);
p0 = mean(n_suc)/50;

[chi2_stat,df,p_value] = sparta_chi_square_test(n_suc,n_tot,p0);


%% Mixtures average time in trap
T_Traps = [T_Traps_Mmix1, T_Traps_Mmix2, T_Traps_Mmix3];
std_Traps = [std_Traps_Mmix1;std_Traps_Mmix2; std_Traps_Mmix3]';


figure;
subplot(1,2,1);
errorbar([0.5,2.5,4.5],T_Traps(1:3,1),std_Traps(1:3,1),'s','linewidth',1,'MarkerSize',8); hold on;
errorbar([1,3,5],T_Traps(1:3,2),std_Traps(1:3,2),'s','linewidth',1,'MarkerSize',8); 
errorbar([1.5,3.5,5.5],T_Traps(1:3,3),std_Traps(1:3,3),'s','linewidth',1,'MarkerSize',8); 
grid on; set(gca,'XTick',[1,3,5],'XTickLabel',ID,'XLim',[0, 6]);
ylim([-1 11]);
axis square;
ylabel('Time (sec)');
pos = get(gca,'Position');
pause(0.1);
legend('Day 1','Day 2','Day 3','Location','northoutside');
pause(0.1);
set(gca,'Position',pos);

n_tot = repmat(10,3,1);
n_suc = T_Traps(:,3);
p0 = mean(n_suc)/10;

[chi2_stat,df,p_value] = sparta_chi_square_test(n_suc,n_tot,p0);
%%


T_Traps = [T_Traps_SMLmix1([2,1,3]),T_Traps_SMLmix2, T_Traps_SMLmix3];
std_Traps = [std_Traps_SMLmix1([2,1,3]);std_Traps_SMLmix2; std_Traps_SMLmix3]';
subplot(1,2,2);
errorbar([0.5,2.5,4.5],T_Traps(1:3,1),std_Traps(1:3,1),'s','linewidth',1,'MarkerSize',8); hold on;
errorbar([1,3,5],T_Traps(1:3,2),std_Traps(1:3,2),'s','linewidth',1,'MarkerSize',8); 
errorbar([1.5,3.5,5.5],T_Traps(1:3,3),std_Traps(1:3,3),'s','linewidth',1,'MarkerSize',8); 
grid on; set(gca,'XTick',[1,3,5],'XTickLabel',{'L';'M';'S'},'XLim',[0, 6]);
ylim([-1 11]);
axis square;

mysetsize(gcf,10,18);

n_tot = repmat(10,3,1);
n_suc = T_Traps(:,3);
p0 = mean(n_suc)/10;

[chi2_stat,df,p_value] = sparta_chi_square_test(n_suc,n_tot,p0)

%% Multiple particle in one spectrum
cd(path_fun);
pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-08\S32P6Mmix1000\short_acquisitions';
[X,~, ~, ~] = sparta_import_shortacq(pat_data); % data with mix of tags
Xp = sparta_preprocess(X,wn);

figure;
plot(wn,Xp(3,:),'LineWidth',1);
not_so_tight('Y');
grid on;

%% particle falling out of trap
cd(path_fun);
pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-08\S31P3Syne0500\short_acquisitions';
[X,~, I, it] = sparta_import_shortacq(pat_data); % data with mix of tags
Xp = sparta_preprocess(X,wn);



% make PCA
LV = 1;
[T,P,expvar] = nipals_pca(mncn(Xp),LV);


figure;
subplot(2,1,1);
mycolor_scatter([261:280]',T(261:280),repmat([1:it]',1,1));
xlabel('Time (sec)')
ylabel('Score LV1 (92.5%)');

hs2 = subplot(2,1,2);
mycolor(wn,Xp(261:280,:),[1:20]');
ylim([0 1300]);
xlabel('Raman shift(cm^-^1)');
ylabel('Intensity');
% pos = get(hs2,'Position');
% pos(2) = pos(2)+0.07;
hcb = colorbar;
set(hcb,'Location','southoutside');
ylabel(hcb, ['Time [sec]'])
caxis([0 it]);
pause(0.1);
% set(hs2,'Position',pos);
mysetsize(gcf,15,20);
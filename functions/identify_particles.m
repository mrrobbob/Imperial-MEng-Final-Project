% clear
close all; clc
pat_fun = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\code\functions';
% data with mix of tags
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S19P6MmixStoc\short_acquisitions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-07\S26P6Mmix1000\short_acquisitions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-08\S32P6Mmix1000\short_acquisitions';

% SML
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S20P8SMLm1000\short_acquisitions';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-07\S27P8SMLm1000\short_acquisitions';
pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-08\S33P8SMLm1000\short_acquisitions';

% corresponding solvent measurements
% pat_data_bck = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S19P4WaterBck\short_acquisitions';
% pat_data_bck = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-07\S26P4WaterBck\short_acquisitions';
% pat_data_bck = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-08\S32P4WaterBck\short_acquisitions';
% SML
% pat_data_bck = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S20P7WaterBck\short_acquisitions';
% pat_data_bck = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-07\S27P7WaterBck\short_acquisitions';
pat_data_bck = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-08\S33P7WaterBck\short_acquisitions';


% data with pure yne-tag
% pat_data_yne = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S18P3Syne0500\short_acquisitions';
% pat_data_yne = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-07\S25P3Syne0500\short_acquisitions';
pat_data_yne = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-08\S31P2Syne1000\short_acquisitions';

% data with pure jon-tag
% pat_data_jon = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-02\S17P8Ljon0250\short_acquisitions';
% pat_data_jon = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-07\S23P6Mjon1000\short_acquisitions';
pat_data_jon = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-06-08\S29P6Mjon1000\short_acquisitions';

mncn = @(x) (x-mean(x));

cd(pat_fun);
load spartaalphawn.mat
% load M_matrix_20230602.mat
% load M_matrix_20230607.mat
load M_matrix_20230608.mat
idx = M(:,strcmp(pat_data(end-31:end-19),ID_M)); % index for successfull trap
idx = logical(idx);

%%
clear ID_M M

% variable range for rep-tag
[~,a_rep] = min(abs(wn-2080));
[~,b_rep] = min(abs(wn-2120));

% variable range for jon- and yne-tag
[~,a_jy] = min(abs(wn-2160));
[~,b_jy] = min(abs(wn-2240));

% %% find rep-tag
[X,~, ~, ~] = sparta_import_shortacq(pat_data); % data with mix of tags
[Xbck,~, ~, ~] = sparta_import_shortacq(pat_data_bck); % corresponding data with pure solvent
[~,I] = max(max(X,[],2)); % remove extreme
idx(I) = false;
% figure;
% subplot(2,1,1); plot(wn,X(idx,:));
% subplot(2,1,2); plot(wn,Xbck); shg
Xp = sparta_preprocess(X,wn);
Xbckp = sparta_preprocess(Xbck,wn);

% plot(wn,Xbck)
%%
% estimate upper limit for solvelt residuals (around rep-tag signal) 
idx_bck = false(size(Xbckp,1),1);
idx_bck(randperm(size(Xbckp,1),ceil(size(Xbckp,1)*0.2))) = true;
mnX = mean(Xbckp(idx_bck,a_rep:b_rep)); % mean of background
E = Xbckp(~idx_bck,a_rep:b_rep)-mnX;
uCI = sparta_upperCI_Qstat(E); % this is the upper limit
clear idx_bck E I

% find rep-tag signal
Ex = Xp(:,a_rep:b_rep)-mnX;
Q = diag(Ex*Ex');
idx_rep = false(size(Xp,1),1);
idx_rep(Q>uCI) = true;
clear Q Ex uCI

idx_rep = idx.*idx_rep;
% %% Find yne and jon tage
% estimate upper limit for solvelt residuals (around yne- and jon-tag signal) 
idx_bck = false(size(Xbckp,1),1);
idx_bck(randperm(size(Xbckp,1),ceil(size(Xbckp,1)*0.2))) = true;
mnX = mean(Xbckp(idx_bck,a_jy:b_jy)); % mean of background
E = Xbckp(~idx_bck,a_jy:b_jy)-mnX;
uCI = sparta_upperCI_Qstat(E); % this is the upper limit
clear idx_bck E I

% find yne and jon-tag signal
Ex = Xp(:,a_jy:b_jy)-mnX;
Q = diag(Ex*Ex');
idx_jy = false(size(Xp,1),1);
idx_jy(Q>uCI) = true;
clear Q Ex uCI

idx_jy = logical(idx.*idx_jy);

% %% read pure jon- and yne-tag data
[Xyne,~, ~, ~] = sparta_import_shortacq(pat_data_yne);
Xp_yne = sparta_preprocess(Xyne,wn);

[Xjon,~, ~, ~] = sparta_import_shortacq(pat_data_jon);
% [Xjon,~]=AutoSpike_Matrix(Xjon,wn,1);
Xp_jon = sparta_preprocess(Xjon,wn);
%
% get latent structure of tag-signals
LV = 1;
[~,p_yne,expVar_yne] = nipals_pca(mncn(Xp_yne(:,a_jy:b_jy)),LV);
[~,p_jon,expVar_jon] = nipals_pca(mncn(Xp_jon(:,a_jy:b_jy)),LV);

% %% find yne
Proj = (eye(109)-p_jon*p_jon');
idx_bck = false(size(Xbckp,1),1);
idx_bck(randperm(size(Xbckp,1),ceil(size(Xbckp,1)*0.2))) = true;
mnX = mean(Xbckp(idx_bck,a_jy:b_jy)*Proj); % mean of background
E = Xbckp(~idx_bck,a_jy:b_jy)*Proj-mnX;
uCI = sparta_upperCI_Qstat(E); % this is the upper limit

Ex = Xp(:,a_jy:b_jy)*Proj-mnX;
Q = diag(Ex*Ex');
idx_y = false(size(Xp,1),1);
idx_y(Q>uCI) = true;
clear Q Ex uCI

idx_y = idx_jy.*idx_y;


% %% find jon
Proj = (eye(109)-p_yne*p_yne');
idx_bck = false(size(Xbckp,1),1);
idx_bck(randperm(size(Xbckp,1),ceil(size(Xbckp,1)*0.2))) = true;
mnX = mean(Xbckp(idx_bck,a_jy:b_jy)*Proj); % mean of background
E = Xbckp(~idx_bck,a_jy:b_jy)*Proj-mnX;
uCI = sparta_upperCI_Qstat(E); % this is the upper limit

Ex = Xp(:,a_jy:b_jy)*Proj-mnX;
Q = diag(Ex*Ex');
idx_j = false(size(Xp,1),1);
idx_j(Q>uCI) = true;
clear Q Ex uCI

idx_j = idx_jy.*idx_j;
% %%
idx_all = [idx_rep,idx_j,idx_y];
N_Traps_SMLmix3 = sparta_count_traps(idx_all);
[T_Traps_SMLmix3, std_Traps_SMLmix3] = sparta_time_traps(idx_all);
%%
n_suc = [sum(idx_j);sum(idx_y);sum(idx_rep)];
n_traps = repmat(1000,3,1);
p0 = mean(n_suc)/1000;

[chi2_stat,df,p_value] = sparta_chi_square_test(n_suc,n_traps,p0);



%%

figure;
subplot(2,1,1);

plot(wn(a_jy:b_jy),Xbckp(:,a_jy:b_jy)*Proj,'y'); hold on;
plot(wn(a_jy:b_jy),Xp(5,a_jy:b_jy),'LineWidth',2);
plot(wn(a_jy:b_jy),Xp(5,a_jy:b_jy)*Proj,'LineWidth',2);


subplot(2,1,2);
plot(wn(a_jy:b_jy),p_jon);
shg

% *Xp(:,a_jy:b_jy)


%% MCR
figure
P = [p_yne,p_jon];
plot(wn(a_jy:b_jy),P); %hold on;
% plot(wn(a_jy:b_jy),sum(P,2),'--');
%%
[copt,sopt,~,~,~,~]=als(Xp(idx_jy,a_jy:b_jy),P);


%% plot signal estimates
figure; 
plot(wn(a_jy:b_jy),sopt(1,:),'b'); hold on;
plot(wn(a_jy:b_jy),sopt(2,:),'r');
plot(wn(a_jy:b_jy),P(:,1),'--b');
plot(wn(a_jy:b_jy),P(:,2),'--r');

line([2200, 2200],ylim); % Jon
line([2214,2214], ylim); % yne

%%
LV = 3;
[~,P,expVar] = nipals_pca(mncn(Xp(idx_jy,a_jy:b_jy)),LV);
figure;
plot(wn(a_jy:b_jy),P);
line([2200, 2200],ylim); % Jon
line([2216,2216], ylim); % yne


%%

plot(wn(a_jy:b_jy),Xp_yne(:,a_jy:b_jy)); shg

%%

Xp_test = Xp(idx_jy,a_jy:b_jy);
for i = 1:size(Xp_test,1)
    plot(wn(a_jy:b_jy),Xp_test(i,:));
    line([2200, 2200],ylim); % Jon
    line([2214,2214], ylim); % yne
    pause
end



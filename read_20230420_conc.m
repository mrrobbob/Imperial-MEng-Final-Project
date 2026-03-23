% read data 20230420
% clear; close all; clc;

mncn = @(x) (x-mean(x));

pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-05-24';
% pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-04-20';
pat_fun = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\code\functions';

cd(pat_data);
% fdir = dir('PLGA*');
fdir = dir('New*');
pat_ext = '\short_acquisitions';

pol = 2;
win = 11;
der = 2;
LV = 1;

cd(pat_fun);
figure;
for i = 1:size(fdir,1)
    pat = [pat_data,'\',fdir(i).name,'\',pat_ext];
    [X,ID, I, it] = sparta_import_shortacq(pat);
    Xsavgol = sparta_SavGol_filter(X,pol,win,der);
    [T,~,expvar] = nipals_pca(mncn(Xsavgol),LV);

%     subplot(4,2,i);
    subplot(3,1,i);
    mycolor_scatter([1:size(T,1)']',T*-1,repmat([1:it]',I,1));
    title(strrep(fdir(i).name,'_',' '));
    ylabel(['Score, LV1 (',num2str(expvar(1,1),'%.2f'),'%)']);
    hcb = colorbar;
    caxis([0, it]);
    ylabel(hcb, ['Time [sec]'])
    not_so_tight, grid on;
end

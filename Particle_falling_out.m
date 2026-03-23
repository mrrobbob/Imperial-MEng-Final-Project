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
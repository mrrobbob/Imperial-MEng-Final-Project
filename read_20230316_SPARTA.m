% Read SPARTA 20230316 - 1 sec short acquisition spectra
clear; close all; clc;
% cd('C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\code');
cd('C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-03-16\test_short_acqu\short_acquisitions');

load spartawn.mat

fdir = dir('*.csv');

mncn = @(x) (x-mean(x)); 

% for i = 1:43
%     M = readmatrix(fdir(i).name);
%     mycolor(wn,M(2:end,2:end)',[1:30]'); shg
%     pause
%     close all
% end

m = 2000;
n = 43;
it = 30;

X = nan(m,n*it);
id = nan(n*it,1);

for i = 1:43
    Mtemp = readmatrix(fdir(i).name);
    X(:,it*(i-1)+1:it*i) = Mtemp(2:end,2:end);
    id(it*(i-1)+1:it*i) = repmat(str2num(fdir(i).name(1:end-4)),it,1);
end
[~,I] = sort(id);
id = id(I);
X = X(:,I);

pol = 2;
win = 55;
der = 2;
Xp = sparta_SavGol_filter(X,pol,win,der);

LV = 1;
[T,P,expVar] = nipals_pca(mncn(Xp'),LV);



figure;
h1 = subplot(3,1,1);
plot(wn,X); not_so_tight('Y'); grid on;
h2 = subplot(3,1,2);
plot(wn,Xp);  not_so_tight('Y'); grid on;
h3 = subplot(3,1,3);
plot(wn,P); not_so_tight('Y'); grid on;
linkaxes([h1,h2,h3],'x')

figure;
% plot(1:size(X,2),T*-1,'.','MarkerSize',10);
mycolor_scatter([1:size(X,2)']',T*-1,repmat([1:30]',n,1));
hcb = colorbar;
caxis([1 30]);
ylabel(hcb, ['Time [Sec]'])
not_so_tight, grid on;

%%
figure
mycolor(wn,X(:,61:90)',[1:30]'); shg
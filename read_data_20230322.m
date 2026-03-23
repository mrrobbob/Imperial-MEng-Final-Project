cd('C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\B-Raman\20230322');
% close all; clear; clc
fdir = dir('*DMS*');
% fdir = dir('*THF*');
J = 2000;
I = size(fdir,1);
X = nan(I,J);
ID = cell(I,1);

for i = 1:I
    Mtemp = readmatrix(fdir(i).name); 
    X(i,:) = Mtemp(3,:);
    ID{i} = fdir(i).name;
end
exwl = 785;
wl = Mtemp(2,:);
wn = 1e7*(1/exwl-1./wl);
clear i Mtemp

figure;
subplot(2,1,1);
plot(wn,X); grid on; not_so_tight('Y');
title('DMSO');
% title('THF');
xlabel('Raman shift (cm^-^1');

subplot(2,1,2);
plot(wl,X); grid on; not_so_tight('Y');
xlabel('Wavelenght (nm)')
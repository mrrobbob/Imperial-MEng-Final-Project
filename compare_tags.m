clear; close all; clc;

pat1 = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-03-29\';
pat2 ='\short_acquisitions';
folder = [{'PLGA_3yne_stock'};...
            {'PLGA_cyano_stock_1'};...
            {'PLGA_EdU_stock'}];

cd('C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\code\functions');
[X_3yne,~,~,~] = sparta_import_shortacq([pat1,folder{1},pat2]);
[X_cyano,~,~,~] = sparta_import_shortacq([pat1,folder{2},pat2]);
[X_EdU,~,~,~] = sparta_import_shortacq([pat1,folder{3},pat2]);
cd(pat1);


load spartawn.mat

figure;
h(1) = subplot(3,1,1); plot(wn,X_3yne); not_so_tight('Y'); grid on;
h(2) = subplot(3,1,2); plot(wn,X_cyano); not_so_tight('Y'); grid on;
h(3) = subplot(3,1,3); plot(wn,X_EdU); not_so_tight('Y'); grid on;
linkaxes([h(1),h(2),h(3)],'x');
shg
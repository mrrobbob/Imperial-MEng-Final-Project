clear; close all; clc;

pat1 = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-04-20\';
pat2 ='\short_acquisitions';
pat3 = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-03-29\PLGA_3yne_stock';
folder = [{'PLGA_NPs_2yne_150mlstock_150mlPBS'};...
            {'PLGA_NPs_2yne_150mlstock_200mlPBS'};...
            {'PLGA_NPs_2yne_150mlstock_350mlPBS'}];

cd('C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\code\functions');
[X_2yne,~,~,~] = sparta_import_shortacq([pat1,folder{2},pat2]);
[X_3yne,~,~,~] = sparta_import_shortacq([pat3,pat2]);
cd(pat1);


load spartaalphawn.mat

figure;
h(1) = subplot(2,1,1); plot(wn,X_2yne); not_so_tight('Y'); grid on;
title('2yne tag');
h(2) = subplot(2,1,2); plot(wn,X_3yne); not_so_tight('Y'); grid on;
title('3yne tag');
linkaxes([h(1),h(2)],'x');
shg
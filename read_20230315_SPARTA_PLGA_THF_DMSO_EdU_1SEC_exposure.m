clear; close all; clc;
cd('C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-03-15\PLGA_THF_DMSO_EdU_1SECEXP');
[NUM, TXT] = xlsread('measurement.csv');

X = NUM(24:end,:);
wn = NUM(2,:);

plot(wn,X);
%%
wn_new = wn;
load spartawn.mat;
wn_new = wn_new';


isequal(wn,wn_new)

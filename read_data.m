clear; close all; clc

idx_date = 11;

path_fun = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\code\functions';
path_main = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman';
cd(path_main)
path_date = dir('2023*');
cd([path_main,'\',path_date(idx_date).name]);
path_exp = dir('S*');
path_acq = 'short_acquisitions';
Bck_str = 'WaterBck';

cd(path_fun);
load spartaalphawn.mat

M = nan(1000,size(path_exp,1)); % collects info whether a particle is present or not
ID_M = cell(size(path_exp,1)); % Experiment ID
for j = 1:size(path_exp,1)
    % check that current (j) experiment is not background
    if strcmp(path_exp(j).name(6:13),Bck_str)
        continue
    end
    ID_M{j} = path_exp(j).name; % ID of experiment
    % find the correct background for current (j) experiment
    for i = 1:size(path_exp,1)
        if strcmp([path_exp(j).name(1:3),Bck_str],path_exp(i).name([1:3,6:13]))
%             read background measurements
            [Xbck,~,~,~] = sparta_import_shortacq([path_main,'\',path_date(idx_date).name,'\',path_exp(i).name,'\',path_acq]);
            Xbck(:,1980) = Xbck(:,1980)+1E4; % add spike
%             subplot(2,1,1); plot(wn, Xbck); title(num2str(i));
            [Xbck,~]=AutoSpike_Matrix(Xbck,wn,1);
%             subplot(2,1,2); plot(wn, Xbck); title(num2str(i));
            Xbckp = sparta_preprocess(Xbck,wn);
            
            % estimate upper limit for background residuals
            idx_calbck = false(size(Xbckp,1),1);
            idx_calbck(randperm(size(Xbckp,1),ceil(size(Xbckp,1)*0.2))) = true;
            mnX = mean(Xbckp(idx_calbck,:)); % mean of background
            E = Xbckp(~idx_calbck,:)-mnX;
            uCI = sparta_upperCI_Qstat(E); % this is the upper limit
            clear i idx_calbck E
            break
        end
    end
    
    % sample measurements
    [X,~,~,~] = sparta_import_shortacq([path_main,'\',path_date(idx_date).name,'\',path_exp(j).name,'\',path_acq]);
    X(:,1980) = X(:,1980)+1E4;
    if j == 3 && idx_date == 11
        X(309,499:506) = mean(X(309,[488:498,507:517]));
    end
    subplot(2,1,1);plot(wn,X); title(num2str(j));
    [X,~]=AutoSpike_Matrix(X,wn,1);
    subplot(2,1,2);plot(wn,X); title(num2str(j));
    Xp = sparta_preprocess(X,wn);
    E = Xp-mnX;
    Q = diag(E*E');
    idx = false(size(Xp,1),1);
    idx(Q>uCI) = true;
    M(:,j) = idx;
end
ID_M = ID_M(~isnan(M(1,:)));
M = M(:,~isnan(M(1,:)));
clear; close all; clc;

mncn = @(x) (x-mean(x));        % anonymous function for mean centering (used for PCA) 

% path for functions
pat_fun = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\code\functions';
% path for data
pat_data = 'C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\data\Raman\2023-05-31\S03P8Ljon1000\short_acquisitions';

cd(pat_fun);
load spartaalphawn.mat
[X,ID,I,it] = sparta_import_shortacq(pat_data);
% I: number traps
% it: numer of measurements pr trap

% Asymmetric least squares fit (for background estimation)
lambda_x = 1e4;     % smoothing parameter (generally 1e5 to 1e8) 
p_x = 0.0001;       % asymmetry parameter (generally 0.001)
d_x = 2;            % order of differences in penalty (generally 2)
Xbck = zeros(size(X));
for i = 1:size(X,1)
    Xbck(i,:) = [asysm(X(i,:)', lambda_x, p_x, d_x)]'; % Asymetric least squares to estimate baseline
%     plot(wn,X(i,:)-Xbck(i,:)); hold on;
%     plot(wn,X(i,:));
%     plot(wn,Xbck(i,:)); hold off;
%     pause
end
%%
% Savitzgy Golay filter (for noise filtering)
pol = 2;
win = 11;
der = 0;    % only noise filtering
Xp = sparta_SavGol_filter(X-Xbck,pol,win,der);

% plot all data
figure;
subplot(2,1,1);
plot(wn, X');
not_so_tight('Y'); grid on;

subplot(2,1,2);
plot(wn,Xp');
not_so_tight('Y'); grid on;
shg


%%
% plot measurements from one trap, color by time
trap = 3;   % trap 3
idx = ID==trap;     % index for trap 3
time = [1:it]';
figure;
mycolor(wn,Xp(idx,:),time);
xlabel('Raman shift [cm^-^1]');
ylabel('Intensity')

% add colorbar
hcbar = colorbar;
caxis([min(time) max(time)]);
ylabel(hcbar,'Time [sec]');
%%
% make PCA
LV = 1;
[T,P,expvar] = nipals_pca(mncn(Xp),LV);


figure;
% mycolor_scatter([1:size(T,1)']',T*-1,repmat([1:it]',I,1));
mycolor_scatter([1:size(T,1)']',T,repmat([1:it]',I,1));

ylabel(['Score, LV1 (',num2str(expvar(1,1),'%.2f'),'%)']);
hcb = colorbar;
caxis([min(time) max(time)]);
ylabel(hcb, ['Time [sec]'])
not_so_tight, grid on;


% adjust size of the figure
hight = 12; % cm
width = 16; % cm
mysetsize(gcf,hight,width)
shg

%%
t_start = 1141;
t_end = t_start+it-1;


figure;
ha = axes; grid on; box on; hold on;
for i = t_start:t_end
    plot(ha,wn,Xp(i,:)');
    title(ha,num2str(i));
    pause;
end


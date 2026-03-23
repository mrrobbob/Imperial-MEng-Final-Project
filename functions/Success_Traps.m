clear; close all; clc
cd('C:\Users\ceskilds\OneDrive - Imperial College London\Documents\mydoc20220808\projects\MSc\Chris\code\functions');
load M_matrix_20230602.mat

N_Traps = nan(50,size(M,2)); % number of successful traps

for j = 1:size(M,2)
    count = 0;
for i = 1:20:size(M,1)
    count = count+1;
    if sum(M(i:i+19,j)) > 0
        N_Traps(count,j) = 1;
    else
        N_Traps(count,j) = 0;
    end
end
end
clear i j count

sum(N_Traps)

% This is just a draft
T_Traps = nan(50,size(M,2)); % time in the trap
for j = 1:size(M,2)
    count = 0;
    for i = 1:20:size(M,1)
        count = count+1;
        for k = 0:19
            if M(i+k,j) == 1
                if k<10
                    T_Traps(count,j) = sum(M(i+k+1:i+k+10,j));
                    break
                elseif k>11 && M(i+19,j) == 0
                    T_Traps(count,j) = sum(M(i+k+1:i+19,j));
                    break
                end
            end
        end
    end
end


%%

[phat,CI] = sparta_CI_trap_prob(repmat(50,size(N_Traps,2),1),sum(N_Traps)');

%%
figure;
% plot(1,phat(2),'s','linewidth',2);
xlab = cell(size(ID_M,2),1);
for i = 1:size(xlab,1);
    xlab{i} = [ID_M{i}(6),ID_M{i}(10:13)];
end

errorbar(phat,phat-CI(:,1),'s','linewidth',1,'MarkerSize',10); grid on;
set(gca,'XTick',[1:size(xlab,1)],'XTickLabel',xlab);
ylabel('probability');

shg

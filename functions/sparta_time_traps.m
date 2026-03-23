function [T_Traps, std_Traps] = sparta_time_traps(M)


% This is just a draft
Traps = nan(50,size(M,2)); % time in the trap
for j = 1:size(M,2)
    count = 0;
    for i = 1:20:size(M,1)
        count = count+1;
        for k = 0:19
            if M(i+k,j) == 1
                if k<10
                    Traps(count,j) = sum(M(i+k+1:i+k+10,j));
                    break
                elseif k>11 && M(i+19,j) == 0
                    Traps(count,j) = sum(M(i+k+1:i+19,j));
                    break
                end
            end
        end
    end
end
T_Traps = nan(size(M,2),1);
for i = 1:size(T_Traps,1)
    T_Traps(i) = mean(Traps(~isnan(Traps(:,i)),i));
    std_Traps(i) = std(Traps(~isnan(Traps(:,i)),i));
end

for i = 1:size(T_Traps,1)
    if T_Traps(i)+std_Traps(i)>10
        std_Traps(i) = 10-T_Traps(i);
    end
    if T_Traps(i)-std_Traps(i)<0
        std_Traps(i) = T_Traps(i);
    end
end


end
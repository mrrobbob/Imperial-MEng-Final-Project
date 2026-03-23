function N_Traps = sparta_count_traps(M)

Traps = nan(50,size(M,2)); % number of successful traps

for j = 1:size(M,2)
    count = 0;
for i = 1:20:size(M,1)
    count = count+1;
    if sum(M(i:i+19,j)) > 0
        Traps(count,j) = 1;
    else
        Traps(count,j) = 0;
    end
end
end
clear i j count

N_Traps = sum(Traps);
end
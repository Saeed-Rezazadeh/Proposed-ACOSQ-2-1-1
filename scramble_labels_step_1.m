function [codebook , T_1] = scramble_labels_step_1(f , Pr , T_1 , numLevel , bit_pattern)
% Scramble the partition indexes. 
j = 1 ;
remaining_bit = zeros(log2(numLevel) - length(bit_pattern) , 1) ;
for i = 1 : log2(numLevel)
    if (ismember(bit_pattern , i) == zeros(1 , length(bit_pattern)))
        remaining_bit(j) = i ;
        j = j + 1 ;
    end
end

for u_index = 1 : length(T_1)
    primary_x = T_1(u_index , 2) ;
    binary_primary_x = de2bi(primary_x - 1 , log2(numLevel) , 'left-msb') ;
    
    
    binary_secondary_x = binary_primary_x([bit_pattern remaining_bit']) ;
    secondary_x = bi2de(binary_secondary_x , 'left-msb') + 1;
    T_1(u_index , 2) = secondary_x ;
end

% Find the codebook based on the new partition indexes. Note that we only
% change/scramble the partition indexes, therefore, only the order of
% codewords in the codebook changes and the values of the codewords are not modified. 

parfor y = 1 : numLevel
    numerator = 0 ;
    denominator = 0 ;
    for x = 1 : numLevel
        u_index = find (T_1 (: , 2) == x) ;
        u = T_1(u_index , 1) ;
        
        numerator = numerator + Pr (x , y) * sum (u .* f(u_index)) ;
        denominator = denominator + Pr (x , y) * sum (f(u_index)) ;
    end
    codebook(y) = numerator / denominator ;
end
end
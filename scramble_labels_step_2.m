function [codebook , T_2] = scramble_labels_step_2(f , Pr , T_2 , numLevel , bit_index)
for u_index = 1 : length(T_2)
    for y_1_2 = 1 : 4
        primary_x = T_2(u_index , 1 + y_1_2 ) ;
        binary_x = de2bi(primary_x - 1 , log2(numLevel) , 'left-msb') ;
        
        hold_var = binary_x(bit_index(y_1_2)) ;
        binary_x(bit_index(y_1_2)) = binary_x(3) ;
        binary_x(3) = hold_var ;
        
        secondary_x = bi2de(binary_x , 'left-msb') + 1;
        T_2(u_index , 1 + y_1_2) = secondary_x ;
    end
end

codebook = zeros(numLevel , 1) ; 
for y_1_2 = 1 : 4
    for y_prime = 1 : 4
        y = (y_1_2 - 1) * 4 + y_prime ;
        numerator = 0 ;
        denominator = 0 ;
        for x = 1 : numLevel
            u_index = find (T_2 (: , 1 + y_1_2) == x) ;
            u = T_2(u_index , 1) ;
            
            numerator = numerator + Pr (x , y) * sum (u .* f(u_index)) ;
            denominator = denominator + Pr (x , y) * sum (f(u_index)) ;
        end
        codebook(y) = numerator / denominator ;
    end
end
end
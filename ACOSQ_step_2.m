function [SDR_2 , T_2 , codebook] = ACOSQ_step_2(f , Pr , Pr_z , codebook , T_1 , numLevel , delta)
FileID = fopen ('Results.txt' , 'a') ;
D_2 = [2 1] ;
Threshold = 0.001 ;
while ((D_2(1) - D_2(2)) / D_2(2)) > Threshold /4
    D_2(1) = D_2(2) ;
    %% Optimal partitions
    T_u = zeros(length(T_1) , 4) ;
    for u_index = 1 : length(T_1)
        summation  = 0 ;
        d_2 = zeros(4 , 4) ;
        
        u = T_1(u_index , 1) ;
        
        hold_x = T_1(u_index , 2) ;
        binary_x = de2bi(hold_x - 1 , log2(numLevel) , 'left-msb') ;
        x_1 = binary_x(1) + 1 ;
        x_2 = binary_x(2) + 1 ;
        x_1_2 = (x_1 - 1) * 2 + x_2 ;
        
        for y_1 = 1 : 2
            for y_2 = 1 : 2
                y_1_2 = (y_1 - 1) * 2 + y_2 ;
                for x_3 = 1 : 2
                    for x_4 = 1 : 2
                        x_3_4 = (x_3 - 1) * 2 + x_4 ;
                        for y_3 = 1 : 2
                            for y_4 = 1 : 2
                                y = (y_1_2 - 1) * 4 + (y_3 - 1) * 2 + y_4 ;
                                summation = summation + ...
                                    Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1)...
                                    * Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1) ...
                                    * (u - codebook(y)) .^ 2 ;
                            end
                        end
                        d_2(y_1_2 , x_3_4) = summation ;
                        summation = 0 ;
                    end
                end
            end
        end
        [~ , partition_index] = min(d_2 , [] , 2) ;
        T_u (u_index , :) = (x_1_2 - 1) * 4  + partition_index ;
    end
    T_2 = cat(2 , T_1(: , 1) , T_u) ;
    
    %% Optimal centroids
    for y_1_2 = 1 : 4
        for y_3_4 = 1 : 4
            y = (y_1_2 - 1) * 4 + y_3_4 ;
            numerator = 0 ;
            denominator = 0 ;
            
            for x = 1 : numLevel
                u_index = find (T_2(: , 1 + y_1_2) == x) ;
                u = T_2(u_index , 1 ) ;
                
                numerator = numerator ...
                    + Pr(x , y) * sum(u .* f(u_index)) ;
                
                denominator = denominator ...
                    + Pr(x , y) * sum(f(u_index)) ;
            end
            codebook(y) = numerator / denominator ;
        end
    end
    %% Distortion
    [D_2(2)] = distortion_2(f , Pr, T_2 , codebook , numLevel , delta) ;
    
    fprintf (FileID , 'Overall D_2 = %f\n' ,D_2(2)) ;
end
fprintf (FileID , 'Overall D_2 = %f\n' ,D_2(2)) ;
SDR_2 = 10 * log10(1 / D_2(2)) ;
fprintf (FileID , 'Overall SDR_2 = %f\n' , SDR_2 ) ;

fprintf (FileID , '=================\n') ;
fclose (FileID) ;

end
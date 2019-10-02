function [SDR , Distortion , T , codebook] = ACOSQ_step_1 (f , Pr , numLevel , T , codebook , delta , b)

FileID = fopen ('Results.txt' , 'a') ;

D = [1 2] ;

while  abs ((D(2) - D(1)) / D(2)) >= (0.001 /4)
    D(1) = D(2) ;
    %% Optimal partitions
    parfor u_index = 1 : length(T)
        summation = 0 ;
        d_1 = zeros (1 , numLevel) ;
        u = T(u_index , 1) ;
        for i = 1 : numLevel
            for j = 1 : numLevel
                summation = summation + Pr(j , b(i)) * (u - codebook (j)) ^ 2 ;
            end
            d_1 (i) = summation ;
            summation = 0 ;
        end
        [~ , partition_index] = min (d_1) ;
        T_u (u_index) = partition_index ;
    end
    T (: , 2) = T_u ;
    %% Optimal Centroids
    parfor j = 1 : numLevel
        numerator = 0 ;
        denominator = 0 ;
        for i = 1 : numLevel
            u_index = find (T (: , 2) == i) ;
            u = T(u_index , 1) ;
            
            numerator = numerator + Pr (j , b(i)) * sum (u .* f(u_index)) ;
            denominator = denominator + Pr (j , b(i)) * sum (f(u_index)) ;
        end
        codebook(j) = numerator / denominator ;
    end
    
    D(2) = distortion_1 (f , numLevel , codebook , delta , Pr , T , b) ;
    fprintf (FileID , 'Overall D_1 = %f\n' ,D(2)) ;
end
SDR = 10 * log10(1 / D (2)) ;
Distortion = D(2);
fprintf (FileID , 'SDR_1 = %f\n' , SDR) ;
fprintf (FileID , '=================\n') ;
fclose (FileID) ;
end
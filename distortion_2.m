function [D_2] = distortion_2(f , Pr , T_2 , codebook , numLevel , delta)
%% Overall adaptive distortion at step 2
summation = 0 ;
parfor x = 1 : numLevel
    for y_1_2 = 1 : 4
        for y_3_4 = 1 : 4
            y = (y_1_2 - 1) * 4 + y_3_4 ;
            
            u_index = find (T_2(: , 1 + y_1_2) == x);
            u = T_2(u_index , 1) ;
            
            summation = summation + delta * Pr(x , y) * sum(f(u_index) .* (u - codebook(y)).^2) ;
            
        end
    end
end
D_2 = summation ;
end
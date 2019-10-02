function [D_4] = distortion_3 (f , Pr , T_3 , codebook , numLevel , delta)

%% Overall adaptive distortion at step 4
summation = 0 ;
parfor x = 1 : numLevel
    for y_1_y_2_y_3 = 1 : 8
        u_index = find (T_3(: , 1 + y_1_y_2_y_3) == x ) ;
        u = T_3(u_index , 1) ;
        for y_4 = 1 : 2
            y = (y_1_y_2_y_3 - 1) * 2 + y_4 ;

            summation = summation  + Pr(x , y) * delta * sum(f(u_index) .*(u - codebook(y)).^2 );
        end
    end
end
D_4 = summation ;
end
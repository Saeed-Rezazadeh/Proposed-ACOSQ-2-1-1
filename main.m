%% The script corresponds to the Algorithm 5 with r = ( 2 1 1)

clc;
clear ;
close all ;

%% The Results.txt
% This file contains the ultimate SDR values for different channel parameters delta and epsilon.
% Also, since the proposed ACOSQ is an iterative algorithm, the distirtion
% value for every iteration is provided in this file for given a epsilon and delta.
FileID = fopen ('Results.txt' , 'a') ;
alpha = 500000 ;



%% Number of Quantization level
numLevel = 16 ;

%% Channel's cross-over probability epsilon
epsilon = unique ([10 ^ -6 10^-5 : 2 * 10^-5 : 10^-4 , 10^-4  : 10^-4  : 10^-3 , 10^-3 , 0.005 0.01 0.05 0.1]);

% Since the design converges to a locally optimal solution, to avoid
% bad local optimums, we use a increase-decrease method.
SIZE = length(epsilon) ;
noise =  [1 : SIZE , SIZE : -1 : 1 , 1 : SIZE , SIZE : -1 : 1] ;


%% SDR parameters
SDR_3 = zeros(length(noise) , 1) ;
SDR_2 = zeros(length(noise) , 1) ;
SDR_1 = zeros(length(noise) , 1) ;
Final_SDR_2 = zeros(12 , length(noise)) ;
Final_SDR_3 = zeros(192 , length(noise)) ;
%% Noise correlation
% The variable delta determines the amount of noise correlation.
for delta = [0 5 10]
    % Set the channel's cross-over probability
    for k = 1 : length(noise)
        i = noise(k);
        
        Pr_1 = [1 - epsilon(i) , epsilon(i) ;
            epsilon(i) , 1 - epsilon(i)] ;
        
        Pr_z = [(1 - epsilon(i) + delta) / (1 + delta)  , epsilon(i) / (1 + delta) ;
            (1 - epsilon(i)) / (1 + delta)  , (epsilon(i) + delta) / (1 + delta)] ;
        
        % Find the channels transition distribution for a given number of
        % quantization levels numLevel
        Pr = Channel_with_Memory(numLevel , epsilon(i) , delta) ;
        
        % Set up the parameters for computing the integrals. Note that in
        % this script all integrals are computed numerically using the
        % Riemann summation.
        delta_u = 8 / 2^ 11 ;
        T_1(: , 1) = -4 : delta_u : 4 ;
        u = T_1(: , 1) ;
        % Compute the source pdf. We herein consider a zero-mean
        % unit-variance Gaussian source distribution.
        f =  1 ./ (sqrt (2 .* pi)) .* exp (-u .^ 2 ./ 2) ;
        f = f ./ (sum(f) .* delta_u) ;
        counter_1 = 0 ;
        counter_2 = 0 ;
        
        % As noted in the Thesis, the ultimate codebook otabined in the
        % last step of the ACOSQ decribed in Section 4.1 is used as the
        % initial state of the proposed ACOSQ with identical bit allocation, noise
        % correlation, and smallest cross-over probability.
        if (k == 1)
            LOAD = ['ACOSQ_2_1_1_delta_' num2str(delta)] ;
            load (LOAD) ;
            codebook_1 = hold_codebook ;
        else
            % We slightly increase the channel's cross-over probability,
            % setting the codebook from the system with small epsilon as
            % the initial state of the system with new epsilon.
            load codebook_1
        end
        % The first step the proposed ACOSQ. Design a 4 bit COSQ
        [SDR_1(k) , ~ , T_1 , codebook_1] = ACOSQ_step_1(f , Pr , numLevel , T_1 , codebook_1 , delta_u , 1 : 16) ;
        
        % save the codebook to initialize the system with the next value of
        % epsilon.
        save('codebook_1'  , 'codebook_1') ;
        
        % save the partition set and codebook for the computing the experimental results.
        Data = ['T\T_1_k_' num2str(k) '_delta_' num2str(delta)]  ;
        save(Data , 'T_1' , 'codebook_1') ;
        
        % Exhaustively search for the two bits to transmit over the channel.
        Bit_index_2 = [3 3 3 3 ; 3 3 3 4 ; 3 3 4 3 ; 3 3 4 4 ; 3 4 3 3 ; 3 4 3 4 ; 3 4 4 3 ; 3 4 4 4 ; 4 3 3 3 ; 4 3 3 4 ; 4 3 4 3 ; 4 3 4 4 ; 4 4 3 3 ; 4 4 3 4 ; 4 4 4 3 ; 4 4 4 4]' ;
        for bit_pattern_1 = [1 2 ; 1 3 ; 1 4 ; 2 1 ; 2 3 ; 2 4 ; 3 1 ; 3 2 ; 3 4 ; 4 1 ; 4 2 ; 4 3]'
            counter_1 = counter_1 + 1 ;
            % Scramble the partition indexes based on the two bits chosen for
            % tranmission such that the bits chosen for transmission are
            % always located in the first two places in the channel input sequence.
            [codebook_2 , init_T_1] = scramble_labels_step_1(f , Pr , T_1 , numLevel , bit_pattern_1') ;
            
            % The second step of the proposed ACOSQ design. In this step,
            % the remaining bits i.e., bit 3 4 are generated adaptive to
            % the received bit y_1 y_2 corresponding to the previous
            % transmissions in step one.
            [SDR_2(k) , T_2 , codebook_2] = ACOSQ_step_2(f , Pr , Pr_z , codebook_2 , init_T_1 , numLevel , delta_u) ;
            Final_SDR_2(counter_1 , k) = SDR_2(k) ;
            
            % save the partition set and codebook for the computing the experimental results.
            Data = ['T\T_2_k_' num2str(k) '_delta_' num2str(delta) , '_counter_1_' num2str(counter_1)]  ;
            save(Data , 'T_2' , 'codebook_2') ;
            
            % As mentioned earlier, this script implements the proposed
            % ACOSQ with r = (2 1 1). In the second step, a 4 bit COSQ
            % designed such that the remaining the 2 bits are generated
            % adaptive to the y_1y_2 received over the feedback link where
            % y_1y_2 are the channel outputs corresponding to the two-tuple
            % transmitted in the first step. We herein for every value of y_1y_2 exhaustively search
            % for the best 1 bits out of the generated 2-tuple for transmission.
            for index_2 = 1 : 16
                bit_pattern_2 = Bit_index_2(: , index_2) ;
                counter_2 = counter_2 + 1 ;
                
                % Scramble the quantization cell indexes such that the
                % 1-bit chosen for transmission are alwayes located in
                % the third place in the channel input
                % sequence.
                [codebook_3 , init_T_2] = scramble_labels_step_2(f , Pr , T_2 , numLevel , bit_pattern_2') ;
                [SDR_3(k) , T_3 , codebook_3] = ACOSQ_step_3(f , Pr , Pr_z , codebook_3 , init_T_2 , numLevel , delta_u) ;
                
                
                Final_SDR_3(counter_2 , k) = SDR_3(k) ;
                Data = ['T\T_3_k_' num2str(k) '_delta_' num2str(delta) , '_counter_2_' num2str(counter_2)]  ;
                save(Data , 'T_3' , 'codebook_3') ;
            end
        end
    end
    myFinal_SDR_2 = max(Final_SDR_2 , [] , 1) ;
    final_SDR_2 = zeros(SIZE , 1) ;
    
    % Pick the best SDR value after the end of the so-called
    % increase-decrease method.
    
    % Variable final_SDR_2 provides the SDR values for proposed ACOSQ 2-2
    for i = 1 : SIZE
        index = find (noise == i) ;
        hold_var  = myFinal_SDR_2(index) ;
        hold_var = hold_var (:) ;
        [final_SDR_2(i)  , index] = max(hold_var) ;
        fprintf (FileID , 'i = %d\n' , i) ;
        fprintf (FileID , '\nfinal_SDR_2 = %f' , final_SDR_2(i)) ;
        fprintf (FileID , '\nindex = %d\n' , index) ;
    end
    myFinal_SDR_3 = max(Final_SDR_3 , [] , 1) ;
    final_SDR_3 = zeros(SIZE , 1) ;
    % Variable final_SDR_3 provides the SDR values for proposed ACOSQ 2-1-1
    for i = 1 : SIZE
        index = find (noise == i) ;
        hold_var  = myFinal_SDR_3(index) ;
        hold_var = hold_var (:) ;
        [final_SDR_3(i)  , index] = max(hold_var) ;
        fprintf (FileID , 'i = %d\n' , i) ;
        fprintf (FileID , '\nfinal_SDR_3 = %f' , final_SDR_3(i)) ;
        fprintf (FileID , '\nindex = %d\n' , index) ;
    end
    Data = ['Proposed_ACOSQ_2_1_1_delta_' num2str(delta)] ;
    save (Data , 'final_SDR_3' , 'myFinal_SDR_3' , 'Final_SDR_3' , 'final_SDR_2' , 'myFinal_SDR_2' , 'Final_SDR_2' , 'epsilon') ;
end

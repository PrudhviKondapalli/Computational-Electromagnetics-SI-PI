%% ASSIGNMENT - 1

% Author: PRUDHVI KONDAPALLI
% Email: Prudhvi.kondapalli@colorado.edu
% Credits : Prof. Mohammed Hadi
% Course: ECEN 5524 Principles of Computational Electromagnetics for SI&PI

%% Question - 1
clc;
clear all
%Create a growing matrix array called I
I = [];                                                 
%Declare M=1 intially
M = 1;                                                  
%p_fin is the variable where we store the value from the richardsons
% extrapolation formulae
probable_final_convergent_value = (2/M)*rich_exp(M);    
%Store that p_fin value for M=1 in the matrix 
I(1,1) = probable_final_convergent_value;               
%Increment M by 1, now becomes 2
M = 2;                                                  
%Calculate p_fin for M=2
probable_final_convergent_value = (2/M)*rich_exp(M);    
%Store that p_fin for M = 2 in the same colomn next row 
I(2,1) = probable_final_convergent_value;                        
%Calculate the difference between the two computed values for M=2 and M=1
diff = I(2,1) - I(1,1);                                 
%if difference was lesser than the required accuracy, end the program
%Increment the row by 2
row_val = 2;                                            

                                                        
%Start an while loop untill diff is > 0.0002
while(diff > 0.0002)
    %Multiply M in the order of 2
    M = M*2;                                            
    %Increment row value by 1 so that we place the latest calculated 
    % element in the next row
    row_val = row_val+1;                                
    %Place the new p_fin value with M to the order 2 from previous in the 
    % next cell in the matirx
    I(row_val,1) = (2/M)*rich_exp(M);                   
    %Check for convergence again, and loop back if convergence is not 
    % acheived
    diff = I(row_val,1) - I(row_val-1,1);              
end

%Display M to see where we acheived convergence
disp(M/2)                                               
    

%% Question - 2

clear all
clc
%Create a growing matrix array called C
C = [];                                                 
%Declare a variable count_M, this keeps track of how many times M was been 
% multiplied 
count_M = 0;                                            
%Intially M is 1
M = 1;                                                  
%Calculate p_fin with M=1
probable_final_convergent_value = (2/M)*rich_exp(M);    
%Place the value with M=1 in element 1(colomn 1 and row 1)
C(1,1) = probable_final_convergent_value;                 
%Increment the count_M by 1 as we are gonna multiply M
count_M = count_M + 1;                                  
%Assign that count_M to dec_count(dec_count is decerment counter to apply
dec_count = count_M;                                    
%M becomes 2
M = 2;                                                  
%Recalculate convergence value with M=2        
probable_final_convergent_value = (2/M)*rich_exp(M);    
%Place that in the array matrix in the next row
C(2,1) = probable_final_convergent_value;                  
%Check for convergence acheived
diff = C(2,1) - C(1,1);                                 
%Initialise diagonal_diff to 1
diagonal_diff = 1;                                      
%Now row becomes 2 and col becomes 1 as thats what we are   
row = 2;                                                
col = 1;

%Run a while loop untill convergence achieved
while((diff > 0.0002) && (diagonal_diff > 0.0002))  
    %We want to decrement the counter to apply richardsons formula
    while (dec_count > 0)                               
        %Increment the column value to allow the next new_element 
        % calculation
        col = col+1;                                    
        %Pass the Col, and C array elements required to calulate the 
        % extrapolated value
        C(row,col) = new_element(col, C(row,col-1), C(row-1,col-1));    
        %Diagonal Diff now checks for convergence with latest values
        diagonal_diff = C(row,col) - C(row-1, col-1);  
        %decrement the counter and start again untill dec_count reaches 0
        dec_count = dec_count - 1;
        %If covergence is acheived
        if (diagonal_diff < 0.0002) && (dec_count == 0)                  
            break;                                        
        end
                                
    end

    %If the above loop did not exit without convergence being acheived
    if (diagonal_diff > 0.0002)                       
            %Mutiply M by 2
            M = M*2;                                    
            %Incrment count_M for the next loop
            count_M = count_M + 1;                      
            %Increment row
            row = row+1;                                
            %Reset column to 1, since thats where our base values are
            col = 1;                                    
            %Store new p_fin value    
            C(row,col) = (2/M)*rich_exp(M);             
            %Check for convergnce
            diff = C(row,1) - C(row-1,1);               
            %reset dec_count with incremented M counter value
            dec_count = count_M;                             
    end
end

%% Additional Question

clear all
clc
%Increment in the steps of 0.0001 from 0 to 1
x = 0 : 0.001 : 1;
%Assign all the expressions to the functions
f_1 = 1./sqrt(1-(x.^2));
f_2 = 1./sqrt(1-(x.^2)) - 1./sqrt(2*(1-x));
f_3 = 1./sqrt(2*(1-x));
%Plot the experssions and annotate them properly
plot(x,f_1, 'o',LineWidth=0.5);
hold on
plot(x, f_2, '-.',LineWidth=2);
hold on
plot(x, f_3, '--', LineWidth=2);
hold on
xline(1, '.', Linewidth=2)
xlim([0 1])
ylim([0 5])
ylabel("Function of x [f(x)]")
xlabel("X values")
title("Before and After Singularity Removal")
legend('F_1(x) - Integrand with singularity','F_2(x) - Integrand split ', ...
    'F_3(x) - Singularity contained within 2nd Integral');
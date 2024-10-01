function [Total_integral_value] = integral_solver(fun) %Initial function declaration
%Lets say the total_integral value is 0
Total_integral_value = 0;
M = 32;
syms x
%Run a nested loop to mimic double summation
for x = 1:M
    for x = 1:M
        %calculate the integral value at that moment
        integral_value_at_M = fun;
        %Add it to the total as we need summation
        Total_integral_value = Total_integral_value + integral_value_at_M;
    end
end
end
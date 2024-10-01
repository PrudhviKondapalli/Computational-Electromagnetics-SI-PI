function [Total_integral_value] = rich_exp(M) %Initial function declaration
%Lets say the total_integral value is 0
Total_integral_value = 0;
%Run a nested loop to mimic double summation
for i = 1:M
    for j = 1:M
        %calculate the integral value at that moment
        integral_value_at_M = 1/sqrt(((i-0.5)^2) + ((j-0.5)^2));
        %Add it to the total as we need summation
        Total_integral_value = Total_integral_value + integral_value_at_M;
    end
end
end
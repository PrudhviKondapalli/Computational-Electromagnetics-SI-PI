%% ASSIGNMENT - 2

% Author: PRUDHVI KONDAPALLI
% Email: Prudhvi.kondapalli@colorado.edu
% Credits : Prof. Mohammed Hadi
% Course: ECEN 5524 Principles of Computational Electromagnetics for SI&PI

%% Additional Question - 1
clc;
clear all

% Get the number of iterations from the user
N = input('Enter the value of Iterations I: ');

% Loop through each iteration
for I = 1:N
    % Create a growing matrix array called A
    A = zeros(I, I);
    B = zeros(I, 1);

    syms x

    % Calculate entries of matrix A using symbolic integration
for n = 1:I
    for n = 1:I
            A(n,n) = int(un_x(n), [0 1]);
    end
end


    % Calculate entries of vector B using symbolic integration
    for n = 1:I
        B(n, 1) = -(int(u0_x, [0 1]));
    end

    % Symbolic variable for coefficients C
    syms 'C%d' [N 1]
    % Solve for C directly as a matrix
    C = A \ B;

    % Display matrices A, B, and vector C
    disp('Matrix A:');
    disp(A);

    disp('Vector B:');
    disp(B);

    disp('Matrix C:');
    disp(C);

    % Define symbolic expressions for the first term and second term
    first_u0 = u0_x;
    second_un = 0;

    % Accumulate the sum for all iterations
    for z = 1:I
        second_un = second_un + C(z, 1) * un_x(z);
    end

    % Calculate the Total expression
    Total = first_u0 + second_un;

    % Display the Total expression
    disp('The Total expression is:');
    disp(Total);
    disp(vpa(Total));

    % Define x values
    x_values = [0, 0.2, 0.4, 0.6, 0.8, 1];

    % Evaluate the Total expression for each x value
    for idx_x = 1:length(x_values)
        x_val = x_values(idx_x);
        total_val = vpa(subs(Total, x, x_val));
        x_matrix(idx_x, I) = total_val;
    end
end

% Display the matrix of Total expressions for different x values
disp('Matrix of Total expressions for different x values:');
disp(x_matrix);

% Calculate the Least Squares Error (LSE)
LSE = 0;
for i = 1:length(x_values)
    LSE = LSE + (abs(x_matrix(i, N) - x_matrix(i, N - 1)))^2;
end

LSE = LSE / length(x_values);

% Display the Least Squares Error
disp('Least Squares Error (LSE):');
disp(LSE);

%% Addition 
clc
clear all 

x_values = [0, 0.2, 0.4, 0.6, 0.8, 1];
o1_values = [0,3.388,6.00128,7.31072,6.43029,2];
o2_values = [0,2.0722,3.78752,4.78,4.5,2];
o3_values = [0,1.510186,2.72864,3.48736,3.439146,2];

plot(x_values,o1_values,'LineWidth',2);
hold on
plot(x_values,o2_values,'LineWidth',2);
hold on
plot(x_values,o3_values,'LineWidth',2);
title('Comparison of 1^{st},2^{nd} and 3^{rd} order basis function approximations')
xlabel('x intervals')
ylabel('\Phi(x)')
legend('First order u_0(x)=2x','Second order u_0(x)=2x^2','Third order u_0(x)=2x^3')

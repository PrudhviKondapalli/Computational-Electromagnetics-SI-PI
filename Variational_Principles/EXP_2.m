%% ASSIGNMENT - 2

% Author: PRUDHVI KONDAPALLI
% Email: Prudhvi.kondapalli@colorado.edu
% Credits : Prof. Mohammed Hadi
% Course: ECEN 5524 Principles of Computational Electromagnetics for SI&PI

%% Question - 1
clc;
clear all

% Get the number of iterations from the user
N = input('Enter the value of Iterations I: ');

syms x
% Loop through each iteration
for I = 1:N
    % Create a growing matrix array called A
    A = zeros(I, I);
    B = zeros(I, 1);

    % Calculate entries of matrix A using symbolic integration
    for m = 1:I
        for n = 1:I
             A(m, n) = int((un_prime_x(m)*un_prime_x(n)) + (un_x(m)*un_x(n)), [0 1]);
        end
    end
    
    % Calculate entries of vector B using symbolic integration
    for n = 1:I
        B(n, 1) = 0 - (int((u0_prime_x*un_prime_x(n)) + (u0_x*un_x(n)), [0 1]));
    end

    % Symbolic variable for coefficients C
    syms 'C%d' [I 1]

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
    % Display the matrix of Total expressions for different x values
    disp('Matrix of Total expressions for different x values:');
    disp(x_matrix);

    % Calculate the Least Squares Error (LSE)
    LSE = 0;
    if I > 1
        for i = 1:length(x_values)
            LSE = LSE + (abs(x_matrix(i, I) - x_matrix(i, I - 1)))^2;
        end
        LSE = LSE / length(x_values);
        RMSE = sqrt(LSE);
        % Display the Least Squares Error
        disp('Least Squares Error (LSE):');
        disp(LSE);
        disp('Root mean Squared Error (RMSE):');
        disp(RMSE);
        % Check if RMSE is less than 0.00001
        if RMSE < 0.000002
            disp(['Breaking out at I = ' num2str(I)]);
            disp(['Corresponding x values in the matrix at I = ' num2str(I)]);
            disp(x_matrix(:, I));
            disp('RMSE is less than 0.000002. Stopping the code.');
            break; % Exit the loop
        end
    end
end

% ...

% Define markers for each iteration
markers = {'-o', '-s', '-d', '-^', '-v', '-p'};

% Plot all the graphs on the same figure with unique markers
figure('Position', [100, 100, 800, 400]); % Adjust the position and size

hold on;
for I = 1:N
    plot(x_values, x_matrix(:, I), markers{mod(I, numel(markers)) + 1}, 'LineStyle', '--'); % Use dashed lines
end
hold off;
title('Approximation solutions by Rayleigh Ritz Equation');
xlabel('x values');
ylabel('Approximated value');
legend(strtrim(cellstr(num2str((1:N)', 'N = %d'))), 'Location', 'Best');




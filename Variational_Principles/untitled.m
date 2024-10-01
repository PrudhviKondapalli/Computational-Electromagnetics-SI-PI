function [A, C, Total, result_table] = solve_electromagnetics(N, x_values, tol)
    % Set tolerance for numerical integration
    Tol = tol;

    % Define symbolic functions (if not already defined)
    syms x;
    u0_x = x;
    un_x = @(n) sin(n*pi*x);
    u0_prime_x = diff(u0_x);
    un_prime_x = @(n) diff(un_x(n));

    % Create A matrix
    A = zeros(N, N);
    for m = 1:N
        for n = 1:N
            integrand = (un_prime_x(m)*un_prime_x(n)) + (un_x(m)*un_x(n));
            % Perform numerical integration with specified tolerance
            A(m, n) = integral(integrand, [0 1], 'AbsTol', Tol, 'RelTol', Tol);
        end
    end

    % Create B vector
    B = zeros(N, 1);
    for n = 1:N
        integrand = (u0_prime_x*un_prime_x(n)) + (u0_x*un_x(n));
        % Perform numerical integration with tolerance
        B(n, 1) = integral(integrand, [0 1], 'AbsTol', Tol, 'RelTol', Tol);
    end

    % Solve for C directly as a matrix, considering symbolic nature
    C = sym('C', [N, 1]);
    C = A \ B;

    % Compute the Total expression
    first_u0 = u0_x;
    second_un = 0;
    for i = 1:N
        second_un = second_un + C(i, 1) * un_x(i);
    end
    Total = first_u0 + second_un;

    % Evaluate for different x values
    result_table = table('Size', [length(x_values), 2], 'VariableTypes', {'double', 'double'}, 'VariableNames', {'x', 'Total'});
    for i = 1:length(x_values)
        x_val = x_values(i);
        total_val = vpa(subs(Total, x, x_val));
        result_table(i, :) = {x_val, total_val};
    end

    % Return modified functions and desired outputs
    return;
end

% Example usage (adjust N, x_values, and tolerance as needed)
N = 5;
x_values = [0, 0.2, 0.4, 0.6, 0.8, 1];
tol = 1e-6;  % Set desired tolerance
[A, C, Total, result_table] = solve_electromagnetics(N, x_values, tol);

disp('Matrix A:');
disp(A);

disp('Vector B:');
disp(B);

disp('Matrix C:');
disp(C);

disp('Total expression:');
disp(Total);

disp('Result table:');
disp(result_table);


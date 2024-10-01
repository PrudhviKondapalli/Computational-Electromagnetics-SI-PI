clc;
clear all;

% Values of N to solve for
N_values = [1, 2, 3];

figure(1);
hold on;
grid on;
xlabel("N");
ylabel("α");
title("α vs N");

for N = N_values
    Array_A = zeros(N, N); % Initialize Array_A with appropriate dimensions
    Array_B = zeros(N, 1); % Initialize Array_B with appropriate dimensions

    for k = 1:N
        for j = 1:N
            Ajk = @(x) x.*(1-x.^j).*(-k.*(k+1).*x.^(k-1));
            % Accumulate the integral values in Array_A
            Array_A(j, k) = integral(Ajk, 0, 1);
        end
    end

    for j = 1:N
        Bj = @(x) -(x.^3).*(1-x.^j);
        Array_B(j,1) = integral(Bj, 0, 1);    
    end

    % Solve the system of equations
    Array_alpha = Array_A \ Array_B;

    plot(1:N, Array_alpha, 'o-', 'DisplayName', ['N = ', num2str(N)]);
end

legend('show');
hold off;

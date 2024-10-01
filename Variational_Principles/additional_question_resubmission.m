% Define the basis functions
u0 = @(x) 2*x;
un = @(n, x) x.^n.*(1-x);

% Define the source function
g = @(x) -x;

% Numerical integration function
integrate = @(f, a, b, N) trapz(linspace(a, b, N+1), f(linspace(a, b, N+1)));

% Solve for the coefficients
N_values = [1, 2, 3]; % Number of terms in the series expansion
coeffs = cell(1, length(N_values));

for i = 1:length(N_values)
    N = N_values(i);
    coeffs{i} = zeros(1, N+1);
    for n = 0:N
        fun = @(x) un(n, x) .* g(x);
        coeffs{i}(n+1) = 2 * integrate(fun, 0, 1, 100);
    end
end

% Construct the approximate solutions
x = linspace(0, 1, 100);
phi = cell(1, length(N_values));
for i = 1:length(N_values)
    N = N_values(i);
    phi{i} = u0(x);
    for n = 1:N
        phi{i} = phi{i} + coeffs{i}(n+1) * un(n, x);
    end
end

% Plot the approximate solutions
figure;
hold on;
colors = ['r', 'g', 'b'];
line_styles = {'-', '--', ':'};
legend_entries = cell(1, length(N_values));
for i = 1:length(N_values)
    N = N_values(i);
    plot(x, phi{i}, 'Color', colors(i), 'LineWidth', 2, 'LineStyle', line_styles{i}, 'DisplayName', sprintf('N = %d terms', N));
end
hold off;
xlabel('x');
ylabel('\phi(x)');
title('Approximate Solutions for Poisson''s Equation');
legend('Location', 'best');
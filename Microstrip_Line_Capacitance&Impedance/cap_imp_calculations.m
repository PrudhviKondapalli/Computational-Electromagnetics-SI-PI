%% ASSIGNMENT - 4
% Author: PRUDHVI KONDAPALLI
% Email: Prudhvi.kondapalli@colorado.edu
% Credits : Prof. Mohammed Hadi
% Course: ECEN 5524 Principles of Computational Electromagnetics for SI&PI

%% Question - 1
clc;
clear all

% Get the number of iterations from the user
N = input('Enter value of N: ');
Ratio = input('Enter Width/Height Ratio of the microstrip: ');
height_microstrip = 2;
width_microstrip = Ratio * height_microstrip;
E_r = 4;
Eff_E = ((E_r + 1) / 2) + ((E_r - 1) / (2 * sqrt(1 + (12 * height_microstrip / width_microstrip))));
u_0 = 4 * pi * (10^(-7));          % henry per mm
e_0 = 8.854 * (10^(-12));          % farad per mm
Delta_x = width_microstrip / (2 * N);
beta_max = max((height_microstrip / Delta_x), (200 * pi));

% Pre-compute constants outside the loop
ones_matrix = ones(N, 1);
beta_values = (1:N) - 0.5;
beta_x_h = (beta_values * Delta_x) / height_microstrip;

% Compute the numerator and denominator terms
numerator = cos(beta_x_h) * cos(beta_x_h');
denominator = abs(beta_values) .* (1 + (Eff_E * coth(abs(beta_values))));  % Use element-wise multiplication

% Compute the A_matrix using integral approximation
A_matrix = Delta_x * trapz(beta_values, numerator ./ denominator);

% Solve for Alpha_matrix
Alpha_matrix = A_matrix \ ones_matrix;

% Compute microstrip impedance
V = 1;  % Set V to a numeric value (e.g., 1)
Alpha_matrix = Alpha_matrix * ((pi * e_0 * V) / 2);

Alpha_sum = sum(Alpha_matrix);

Capacitance = (2 * Delta_x / V) * Alpha_sum;
V_p = 1 / sqrt(u_0 * e_0 * Eff_E);
Impedance = 1 / (Capacitance * V_p);

disp(['Microstrip Impedance: ', num2str(double(Impedance))]);

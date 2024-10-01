%% ASSIGNMENT - 4
% Author: PRUDHVI KONDAPALLI
% Email: Prudhvi.kondapalli@colorado.edu
% Credits : Prof. Mohammed Hadi
% Course: ECEN 5524 Principles of Computational Electromagnetics for SI&PI

%% Question - 1
clc;
clear all

% Get the number of iterations from the user
N = input('Enter value of N');
Ratio = input('Enter Width/Height Ratio of the microstrip');
height_microstrip = 2;
width_microstrip = Ratio*height_microstrip;
E_r = 4;
Eff_E = ((E_r+1)/2)+((E_r-1)/(2*sqrt(1+(12*height_microstrip/width_microstrip))));
u_0 = 4*pi*(10^(-7));          % henry per mm
e_0 = 8.854*(10^(-12));          % farad per mm
Delta_x = width_microstrip/(2*N);
beta_max = max((height_microstrip/Delta_x),(200*pi));

% Pre-compute constants outside the loop
ones_matrix = ones(N,1);
syms beta V;  % Declare V as a symbolic variable

% Vectorize beta calculations (if possible for your specific problem)
% beta_xn_h = (beta*(1:N) - 0.5)*Delta_x / height_microstrip;
% beta_xk_h = (beta*(1:N) - 0.5)*Delta_x / height_microstrip;

% Use parfor loop for parallel computation in expensive calculations
A_matrix = zeros(N,N);
parfor n = 1:N
    beta_xn = (beta*(n-0.5)*Delta_x)/height_microstrip;
    for k = 1:N
        beta_xk = (beta*(k-0.5)*Delta_x)/height_microstrip;
        numerator(n,k) = cos(beta_xn)*cos(beta_xk);
        denominator(n,k) = abs(beta)*(1 +(Eff_E*coth(abs(beta))));
    end
end

parfor n = 1:N
  for k =1:N
    A_matrix(n,k) = Delta_x*(int((numerator(n,k)/denominator(n,k)),beta,[0 beta_max]));
  end
end

Alpha_matrix = A_matrix \ ones_matrix;
% V cancels out during calculations
Alpha_matrix = Alpha_matrix*((pi*e_0*V)/2);
Alpha_sum = 0;
for k = 1:N
  Alpha_sum = Alpha_sum + Alpha_matrix(k,1);
end
Capacitance = (2*Delta_x/V)*Alpha_sum;
V_p = 1/sqrt(u_0*e_0*Eff_E);
Impedance = 1/(Capacitance*V_p);
disp(['Microstrip Impedance: ', num2str(double(Impedance))]);

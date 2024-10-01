clear; clc;

% Parameters
mu=4*pi*1e-7; eps=8.853e-12; c0=1/sqrt(mu*eps);
f=1e9; lambda=c0/f; omega=2*pi*f;
dx = lambda/20;
dy = lambda/20;      % Cell sizes (λ/20)
dt = dx/((sqrt(2))*c0);
sinw = 0;
sinbetax = 0;
sinbetay = 0;

angles = 1:90;       % Angle of propagation in degrees
phase_error = zeros(size(angles));

beta_initial = omega/c0;
beta_tilde = beta_initial;
beta_tilde_int = zeros(size(angles));

for i = 1:90
    % %beta_tilde_x = beta_tilde*cos(i*pi/180);
    % %beta_tilde_y = beta_tilde*sin(i*pi/180);
    sinw = sin((omega*(dt/2)))/(dt/2);
    % %sinbetax = sin(beta_tilde_x*(dx/2))/(dx/2);
    % %sinbetay = sin(beta_tilde_y*(dy/2))/(dy/2);

    func = @(beta_tilde) ((mu*eps*((sinw)^2))-((sin((beta_tilde*cos(i*pi/180))*(dx/2))/(dx/2))^2) -((sin((beta_tilde*sin(i*pi/180))*(dy/2))/(dy/2))^2));

    beta_tilde = fsolve(func,beta_initial);
    beta_tilde_int(i) = beta_tilde;

    phase_error(i) = (abs(beta_initial - beta_tilde))*lambda;
end


% Plot the numerical phase error
figure;
plot(angles*pi/180, phase_error, 'b-', 'LineWidth', 2);
xlabel('Angle of Propagation (Degrees)');
ylabel('Phase Error (\lambda)');
title('Numerical Phase Error');

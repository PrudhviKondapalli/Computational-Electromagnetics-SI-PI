clear; clc;

% Parameters
mu = 4*pi*1e-7; eps = 8.853e-12; c0 = 1/sqrt(mu*eps);
f = 1e9; lambda = c0/f; omega = 2*pi*f;
dx = lambda/20;
dy = lambda/20; % Cell sizes (λ/20)

% Different dt values
dt1 = dx/((sqrt(2))*c0); % Original dt
dt2 = dx/((sqrt(2))*2*c0); % dt = dx/((sqrt(2))*2*c0)
dt3 = dx/((sqrt(2))*4*c0); % dt = dx/((sqrt(2))*4*c0)
dt4 = dx/((sqrt(2))*10*c0); % dt = dx/((sqrt(2))*10*c0)

sinw = 0;
sinbetax = 0;
sinbetay = 0;
angles = 1:90; % Angle of propagation in degrees
beta_initial = omega/c0;
beta_tilde = beta_initial;

% Plot for different dt values
figure;
hold on;

for i = 1:90
    sinw = sin((omega*(dt1/2)))/(dt1/2);
    func = @(beta_tilde) ((mu*eps*((sinw)^2))-((sin((beta_tilde*cos(i*pi/180))*(dx/2))/(dx/2))^2) -((sin((beta_tilde*sin(i*pi/180))*(dy/2))/(dy/2))^2));
    beta_tilde = fsolve(func, beta_initial);
    beta_tilde_int1(i) = beta_tilde;
    phase_error1(i) = (abs(beta_initial - beta_tilde))*lambda;
end
plot(angles*pi/180, phase_error1, 'b-', 'LineWidth', 2, 'DisplayName', 'dt = dx/((sqrt(2))*c0)');

for i = 1:90
    sinw = sin((omega*(dt2/2)))/(dt2/2);
    func = @(beta_tilde) ((mu*eps*((sinw)^2))-((sin((beta_tilde*cos(i*pi/180))*(dx/2))/(dx/2))^2) -((sin((beta_tilde*sin(i*pi/180))*(dy/2))/(dy/2))^2));
    beta_tilde = fsolve(func, beta_initial);
    beta_tilde_int2(i) = beta_tilde;
    phase_error2(i) = (abs(beta_initial - beta_tilde))*lambda;
end
plot(angles*pi/180, phase_error2, 'r--', 'LineWidth', 2, 'DisplayName', 'dt = dx/((sqrt(2))*2*c0)');

for i = 1:90
    sinw = sin((omega*(dt3/2)))/(dt3/2);
    func = @(beta_tilde) ((mu*eps*((sinw)^2))-((sin((beta_tilde*cos(i*pi/180))*(dx/2))/(dx/2))^2) -((sin((beta_tilde*sin(i*pi/180))*(dy/2))/(dy/2))^2));
    beta_tilde = fsolve(func, beta_initial);
    beta_tilde_int3(i) = beta_tilde;
    phase_error3(i) = (abs(beta_initial - beta_tilde))*lambda;
end
plot(angles*pi/180, phase_error3, 'g-.', 'LineWidth', 2, 'DisplayName', 'dt = dx/((sqrt(2))*4*c0)');

for i = 1:90
    sinw = sin((omega*(dt4/2)))/(dt4/2);
    func = @(beta_tilde) ((mu*eps*((sinw)^2))-((sin((beta_tilde*cos(i*pi/180))*(dx/2))/(dx/2))^2) -((sin((beta_tilde*sin(i*pi/180))*(dy/2))/(dy/2))^2));
    beta_tilde = fsolve(func, beta_initial);
    beta_tilde_int4(i) = beta_tilde;
    phase_error4(i) = (abs(beta_initial - beta_tilde))*lambda;
end
plot(angles*pi/180, phase_error4, 'k:', 'LineWidth', 2, 'DisplayName', 'dt = dx/((sqrt(2))*10*c0)');

xlabel('Angle of Propagation (Degrees)');
ylabel('Phase Error (\lambda)');
title('Numerical Phase Error');
legend('Location', 'best');
hold off;
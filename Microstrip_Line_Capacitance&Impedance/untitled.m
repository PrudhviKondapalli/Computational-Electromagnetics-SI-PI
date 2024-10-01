% Define the data
dielectric_thickness = [6, 20, 50, 100];  % in mils
diff_impedance = [112.2, 126.68, 128.4, 128.669];

% Create a figure
figure;

% Plot the data
plot(dielectric_thickness, diff_impedance, 'o-');

% Set labels and title
xlabel('Dielectric Thickness (mils)');
ylabel('Differential Impedance (ohms)');
title('Differential Impedance vs. Dielectric Thickness');

% Adjust axis limits if needed
% xlim([0, 120]);  % Set x-axis limits
% ylim([100, 140]); % Set y-axis limits

% Add grid lines
grid on;

% Display the plot
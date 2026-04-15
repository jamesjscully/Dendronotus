% Test Plot Script
set(0, 'DefaultFigureRenderer', 'painters');
clear;
clc;

% Generate some test data
x = linspace(0, 2*pi, 100);
y1 = sin(x);
y2 = cos(x);

% Create a figure
figure;
clf;

% Plot the sine wave
subplot(2, 1, 1);
plot(x, y1, 'r', 'LineWidth', 2);
title('Sine Wave');
xlabel('x');
ylabel('sin(x)');
grid on;

% Plot the cosine wave
subplot(2, 1, 2);
plot(x, y2, 'b', 'LineWidth', 2);
title('Cosine Wave');
xlabel('x');
ylabel('cos(x)');
grid on;

% Save the figure as a file
saveas(gcf, 'test_plot_output.png');

disp('Test plot generated and saved as test_plot_output.png');
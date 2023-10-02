% Kristen Park 
% October 2, 2023
% NGG 6050 
% Homework #9
% Parametric Correlation Coefficient

clear
close all 

%% 1. Plot the relationship between Age and Wing Length.
% Define the data
Age = [3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17];
WingLength = [1.4, 1.5, 2.2, 2.4, 3.1, 3.2, 3.2, 3.9, 4.1, 4.7, 4.5, 5.2, 5.0];

% Plot the relationship between Age and Wing Length.
figure;
plot(Age, WingLength, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
title('Relationship Between Age and Wing Length');
xlabel('Age');
ylabel('Wing Length');

% Add a grid to the plot
grid on;

%% 2. Calculate and plot the regression line.

% Fit a linear regression model
coefficients = polyfit(Age, WingLength, 1);

% Calculate the regression line
regression_line = polyval(coefficients, Age);

% Create the scatter plot
figure;
plot(Age, WingLength, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
hold on;

% Plot the regression line
plot(Age, regression_line, 'r', 'LineWidth', 2);

% Add title and labels
title('Regression Line for Age vs. Wing Length');
xlabel('Age');
ylabel('Wing Length');

% Add a grid to the plot
grid on;

%% 3. Can you reject H0: b = 0? 
% Fit a linear regression model
coefficients = polyfit(Age, WingLength, 1);
slope = coefficients(1);

% Calculate the predicted values
predicted = polyval(coefficients, Age);

% Calculate the residuals
residuals = WingLength - predicted;

% Calculate the sum of squares for regression
SSR = sum((predicted - mean(WingLength)).^2);

% Calculate the sum of squares for residuals (error)
SSE = sum(residuals.^2);

% Calculate the total sum of squares
SST = SSR + SSE;

% Calculate the degrees of freedom
n = length(Age);
df_regression = 1;
df_residual = n - 2;

% Calculate the mean square for regression and residuals
MSR = SSR / df_regression;
MSE = SSE / df_residual;

% Calculate the F-statistic
F_statistic = MSR / MSE;

% Set the significance level
alpha = 0.05;

% Perform the F-test and check if H0 can be rejected
p_value = 1 - fcdf(F_statistic, df_regression, df_residual);
fprintf('F-statistic: %.4f\n', F_statistic);
fprintf('P-value: %.4f\n', p_value);

if p_value < alpha
    fprintf('Reject H0: There is a significant linear relationship between Age and Wing Length.\n');
else
    fprintf('Fail to reject H0: There is no significant linear relationship between Age and Wing Length.\n');
end

%% 4. Calculate and plot the confidence intervals on the slope of the regression.

% Calculate the standard error of the slope (SE)
residuals = WingLength - polyval(coefficients, Age);
MSE = sum(residuals.^2) / (length(Age) - 2); % Mean squared error
SE = sqrt(MSE / sum((Age - mean(Age)).^2)); % Standard error of the slope

% Calculate the t-statistic for a 95% confidence interval
alpha = 0.05;
t_critical = tinv(1 - alpha/2, length(Age) - 2);

% Calculate the confidence interval bounds
slope = coefficients(1);
CI_lower = slope - t_critical * SE;
CI_upper = slope + t_critical * SE;

fprintf('\nSlope: %.4f\n', slope);
fprintf('95%% Confidence Interval for Slope: [%.4f, %.4f]\n', CI_lower, CI_upper);

% Create a scatter plot with the regression line
figure;
plot(Age, WingLength, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
hold on;
plot(Age, polyval(coefficients, Age), 'r', 'LineWidth', 2);

% Plot the confidence interval lines
x_fill = [min(Age), max(Age), max(Age), min(Age)];
y_fill = [polyval(coefficients, min(Age)) - CI_lower, polyval(coefficients, max(Age)) - CI_lower, ...
          polyval(coefficients, max(Age)) + CI_upper, polyval(coefficients, min(Age)) + CI_upper];
fill(x_fill, y_fill, 'g', 'FaceAlpha', 0.3);

% Add title and labels
title('Regression Line with 95% Confidence Intervals for Slope');
xlabel('Age');
ylabel('Wing Length');

% Add a legend
legend('Data', 'Regression Line', '95% CI for Slope', 'Location', 'Northwest');

% Add a grid to the plot
grid on;

%% 5. Calculate r^2 (the coefficient of determination)

% Fit a linear regression model
coefficients = polyfit(Age, WingLength, 1);

% Calculate the predicted values
predicted = polyval(coefficients, Age);

% Calculate the mean of the observed Wing Length
mean_WingLength = mean(WingLength);

% Calculate the total sum of squares (TSS)
TSS = sum((WingLength - mean_WingLength).^2);

% Calculate the residual sum of squares (RSS)
RSS = sum((WingLength - predicted).^2);

% Calculate R-squared (coefficient of determination)
R_squared = 1 - (RSS / TSS);

fprintf('\nR-squared (Coefficient of Determination): %.4f\n', R_squared);

%% 6. Calculate Pearson's correlation coefficient
r = corr(Age', WingLength');

fprintf("\nPearson's Correlation Coefficient (r): %.4f\n", r);

%% 7. Add some noise to the data and see how the regression changes.
Age = [3, 4, 5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 2, 3, 4, 10, 12];
WingLength = [1.4, 1.5, 2.2, 2.4, 3.1, 3.2, 3.2, 3.9, 4.1, 4.7, 4.5, 5.2, 5.0, 6, 5, 3, 1, 2];

% Calculate the standard error of the slope (SE)
residuals = WingLength - polyval(coefficients, Age);
MSE = sum(residuals.^2) / (length(Age) - 2); % Mean squared error
SE = sqrt(MSE / sum((Age - mean(Age)).^2)); % Standard error of the slope

% Calculate the t-statistic for a 95% confidence interval
alpha = 0.05;
t_critical = tinv(1 - alpha/2, length(Age) - 2);

% Calculate the confidence interval bounds
slope = coefficients(1);
CI_lower = slope - t_critical * SE;
CI_upper = slope + t_critical * SE;

% Create a scatter plot with the regression line
figure;
plot(Age, WingLength, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
hold on;
plot(Age, polyval(coefficients, Age), 'r', 'LineWidth', 2);

% Plot the confidence interval lines
x_fill = [min(Age), max(Age), max(Age), min(Age)];
y_fill = [polyval(coefficients, min(Age)) - CI_lower, polyval(coefficients, max(Age)) - CI_lower, ...
          polyval(coefficients, max(Age)) + CI_upper, polyval(coefficients, min(Age)) + CI_upper];
fill(x_fill, y_fill, 'g', 'FaceAlpha', 0.3);

% Add title and labels
title('Noisy Data');
xlabel('Age');
ylabel('Wing Length');

% Add a legend
legend('Data', 'Regression Line', '95% CI for Slope', 'Location', 'Northwest');

% Add a grid to the plot
grid on;

r = corr(Age', WingLength');
fprintf("\nPearson's Correlation Coefficient of noisy data (r): %.4f\n", r);





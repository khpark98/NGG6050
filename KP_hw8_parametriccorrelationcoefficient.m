% Kristen Park 
% September 16, 2023
% NGG 6050 
% Homework #8
% Parametric Correlation Coefficient

clear
close all 

%% Parametric Correlation Coefficient
wingLength = [10.4, 10.8, 11.1, 10.2, 10.3, 10.2, 10.7, 10.5, 10.8, 11.2, 10.6, 11.4];
tailLength = [7.4, 7.6, 7.9, 7.2, 7.4, 7.1, 7.4, 7.2, 7.8, 7.7, 7.8, 8.3];

%% 1. Plot X vs Y. Do they look related? 

% Yes, it looks like there's a positive linear relationship between these 
% two variables. 

scatter(wingLength, tailLength);
xlabel('Wing Length');
ylabel('Tail Length');
title('Wing Length vs. Tail Length');
grid on;


%% 2. Calculate rx,y and ry,x, first using the equations above and then using 
% either the Python numpy funciton corrcoef or Matlab's built-in corrcoef. 
% Did you get the same answers?

% Yes, I got the same answers using the equations and corrcoef.
% 0.8704 for both rxy and ryx! 

%equation
mean_wingLength = mean(wingLength);
mean_tailLength = mean(tailLength);
std_wingLength = std(wingLength);
std_tailLength = std(tailLength);

numerator = sum((wingLength - mean_wingLength) .* (tailLength - mean_tailLength));
denominator = sqrt(sum((wingLength - mean_wingLength).^2) * sum((tailLength - mean_tailLength).^2));
r = numerator / denominator;

numerator_yx = sum((tailLength - mean_tailLength) .* (wingLength - mean_wingLength));
denominator_yx = sqrt(sum((tailLength - mean_tailLength).^2) * sum((wingLength - mean_wingLength).^2));
ryx = numerator_yx / denominator_yx;

fprintf('Pearson correlation coefficient r(x, y) using equation: %.4f\n', r);
fprintf('Pearson correlation coefficient r(y, x) using equation: %.4f\n', ryx);

%corrcoef function 
correlation_matrix = corrcoef(wingLength, tailLength);
rxy = correlation_matrix(1, 2);

correlation_matrix = corrcoef(tailLength, wingLength);
ryx = correlation_matrix(1,2);

fprintf('\nPearson correlation coefficient r(x, y) using corrcoef: %.4f\n', rxy);
fprintf('Pearson correlation coefficient r(y, x) using corrcoef: %.4f\n', ryx);

%% 3. What is the standard error of rx,y? The 95% confidence intervals computed from the standard error?

% The standard error of r(x,y) is 0.1557.

n = length(wingLength); %sample size
SE_r = sqrt((1-r^2)/(n-2));
fprintf('\nStandard Error of r(x, y): %.4f\n', SE_r);

% The 95% CI is [0.5923, 0.9632] 

z = 0.5*log((1+r)/(1-r)); %fisher's z transformation of r
sz = sqrt(1/(n-3)); % st dev of z

alpha = 0.05; %for 95% CI
z_critical = norminv(1 - alpha/2);
margin_of_error = z_critical * sz;
lower = z - margin_of_error;
upper = z + margin_of_error;

r_lower = (exp(2 * lower) - 1) / (exp(2 * lower) + 1);

r_upper = (exp(2 * upper) - 1) / (exp(2 * upper) + 1);
fprintf('95%% Confidence Interval for Pearson''s correlation coefficient: [%.4f, %.4f]\n', r_lower, r_upper);

%% 4. Should the value of rxy be considered significant at the p<0.05 level, 
% given a two-tailed test (i.e., we reject if the test statistic is too 
% large on either tail of the null distribution) for H0: rxy = 0? 

% Yes, it is considered significant at a p<0.05 level.

t_statistic = r * sqrt(n - 2) / sqrt(1 - r^2);
df = n - 2;

% Calculate the critical values for a two-tailed test
alpha_half = alpha / 2;
t_critical_lower = tinv(alpha_half, df);
t_critical_upper = -t_critical_lower; % For the upper tail

% Check if the t-statistic is in the rejection region
if t_statistic < t_critical_lower || t_statistic > t_critical_upper
    fprintf('\nThe correlation is statistically significant (p < %.2f).\n', alpha);
else
    fprintf('\nThe correlation is not statistically significant (p >= %.2f).\n', alpha);
end

%% 5. Yale does the exact same study and finds that his correlation value is 0.75. 
% Is this the same as yours? That is, evaluate H0: r = 0.5. 

% Yes, Yale's correlation value is the same as mine. 

z = 0.5*log((1+r)/(1-r)); %fisher's z transformation of r
rs = 0.75;
z_rs = 0.5*log((1+rs)/(1-rs));

lambda_stat = (z-z_rs)/sqrt(1/(n-3));
df = n-2; 

if lambda_stat < t_critical_lower || lambda_stat > t_critical_upper
    fprintf('\nReject the null hypothesis -- Yale''s correlation value is not the same as yours\n');
else
    fprintf('\nFail to reject the null hypothesis -- Yale''s correlation value is the same as yours\n');
end

%% 6. Finally, calculate the statistical power and sample size needed to reject H0:r=0, when r>=0.5
alpha = 0.05;         % Significance level
power = 0.80;         % Desired power level
effect_size = 0.5;   % Effect size (correlation coefficient under H1)
rho_null = 0;        % Null hypothesis correlation coefficient

z_alpha = norminv(1 - alpha, 0, 1);

% Calculate required sample size
n = (z_alpha + norminv(power, 0, 1))^2 / effect_size^2;

fprintf('\nRequired sample size needed to reject H0:r=0, when r>=0.5, at a 0.8 power level: %d\n', ceil(n));



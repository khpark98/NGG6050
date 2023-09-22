% Kristen Park
% 9/21/23
% Homework #7
% Multiple Comparisons

%%
% First, simulate multiple (say, 1000) t-tests comparing two samples with 
% equal means and standard deviations, and save the p-values. Obviously, 
% at p<0.05 we expect that ~5% of the simulations to yield a "statistically 
% significant" result (of rejecting the NULL hypothesis that the samples 
% come from distributions with equal means).

% Set the parameters for the simulation
num_simulations = 1000;
sample_size = 50; % Size of each sample
alpha = 0.05; % Significance level

% Initialize an array to store the p-values
p_values = zeros(num_simulations, 1);

% Loop for running t-tests and saving p-values
for i = 1:num_simulations
    % Generate two random samples from a standard normal distribution
    sample1 = randn(sample_size, 1);
    sample2 = randn(sample_size, 1);
    
    % Perform a two-sample t-test
    [~, p, ~, ~] = ttest2(sample1, sample2);
    
    % Store the p-value
    p_values(i) = p;
end

% Calculate the proportion of p-values less than alpha
proportion_significant = sum(p_values < alpha) / num_simulations;

% Display the results
fprintf('Proportion of "statistically significant" results: %.2f%%\n', proportion_significant * 100);

%% 
% Second, once you have the simulated p-values, apply both methods to 
% address the multiple comparisons problem.

% Bonferroni Correction
% Calculate the Bonferroni-corrected significance level
alpha_bonf = alpha / num_simulations;

% Initialize an array to store the results of the corrected tests
is_significant = zeros(1, num_simulations);

% Apply the Bonferroni correction
for i = 1:num_simulations
    if p_values(i) <= alpha_bonf
        is_significant(i) = 1; % Reject the null hypothesis (significant)
    else
        is_significant(i) = 0; % Fail to reject the null hypothesis (not significant)
    end
end

sigresults = (is_significant == 1);
proportion_significant = sum(sigresults)/num_simulations;

% Display the results
fprintf('Bonferroni-corrected significance level: %.5f\n', alpha_bonf);
fprintf('Proportion of "statistically significant" results after Bonferroni correction: %.2f%%\n', proportion_significant * 100);

% BH Procedure
% Define the desired False Discovery Rate (FDR)
FDR = 0.05;

% Sort p-values in ascending order
sorted_p_values = sort(p_values);

% Initialize an array to store the results of the corrected tests
is_significant = zeros(1, num_simulations);

% Apply the BH correction
for i = 1:num_simulations
    BH_critical_value = (i / num_simulations) * FDR;
    if sorted_p_values(i) <= BH_critical_value
        is_significant(i) = 1; % Reject the null hypothesis (significant)
    else
        is_significant(i) = 0; % Fail to reject the null hypothesis (not significant)
    end
end

sigresults = (is_significant == 1);
proportion_significant = sum(sigresults)/num_simulations;

% Display the results
fprintf('False Discovery Rate (FDR): %.4f\n', FDR);
fprintf('Proportion of "statistically significant" results after BH procedure: %.2f%%\n', proportion_significant * 100);

%%
% Third, set the sample 1 and sample 2 means to be 1 and 2 respectively, 
% and re-run the exercise. What do you notice? What if you make the 
% difference between means even greater?
% Set the parameters for the simulation
num_simulations = 1000;
sample_size = 50; % Size of each sample
alpha = 0.05; % Significance level

% Initialize an array to store the p-values
p_values = zeros(num_simulations, 1);

% Loop for running t-tests and saving p-values
for i = 1:num_simulations
    % Generate two random samples from a standard normal distribution
    sample1 = randn(sample_size, 1) + 1; %mean = 1
    sample2 = randn(sample_size, 1) + 2; %mean = 2
    
    % Perform a two-sample t-test
    [~, p, ~, ~] = ttest2(sample1, sample2);
    
    % Store the p-value
    p_values(i) = p;
end

% Calculate the proportion of p-values less than alpha
proportion_significant = sum(p_values < alpha) / num_simulations;

% Display the results
fprintf('\nProportion of "statistically significant" results after increasing the difference in means: %.2f%%\n', proportion_significant * 100);
% Bonferroni Correction
% Calculate the Bonferroni-corrected significance level
alpha_bonf = alpha / num_simulations;

% Initialize an array to store the results of the corrected tests
is_significant = zeros(1, num_simulations);

% Apply the Bonferroni correction
for i = 1:num_simulations
    if p_values(i) <= alpha_bonf
        is_significant(i) = 1; % Reject the null hypothesis (significant)
    else
        is_significant(i) = 0; % Fail to reject the null hypothesis (not significant)
    end
end

sigresults = (is_significant == 1);
proportion_significant = sum(sigresults)/num_simulations;

% Display the results
fprintf('Bonferroni-corrected significance level: %.5f\n', alpha_bonf);
fprintf('Proportion of "statistically significant" results after increasing the difference in means and Bonferroni correction: %.2f%%\n', proportion_significant * 100);

% BH Procedure
% Define the desired False Discovery Rate (FDR)
FDR = 0.05;

% Sort p-values in ascending order
sorted_p_values = sort(p_values);

% Initialize an array to store the results of the corrected tests
is_significant = zeros(1, num_simulations);

% Apply the BH correction
for i = 1:num_simulations
    BH_critical_value = (i / num_simulations) * FDR;
    if sorted_p_values(i) <= BH_critical_value
        is_significant(i) = 1; % Reject the null hypothesis (significant)
    else
        is_significant(i) = 0; % Fail to reject the null hypothesis (not significant)
    end
end

sigresults = (is_significant == 1);
proportion_significant = sum(sigresults)/num_simulations;

% Display the results
fprintf('False Discovery Rate (FDR): %.4f\n', FDR);
fprintf('Proportion of "statistically significant" results after BH procedure: %.2f%%\n', proportion_significant * 100);

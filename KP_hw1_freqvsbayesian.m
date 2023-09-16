% Kristen Park 
% September 5, 2023
% NGG 6050 
% Homework #1

clear
close all 

%% Exercise 1
%If someone gets a positive test, is it "statistically significant" at the p<0.05 level? 
%Why or why not

%In this situation, it would be inappropriate to use the p-value. The
%p-value is the probability of obtaining the results R given the null 
%hypothesis H0. It is generally used in the context of
%statistical hypothesis testing. In this case, there are clear false positive
%and false negative rates that you can use to calculate the number of
%expected results, so hypothesis testing is not necessary to understand the
%performance of the test. 

%% Exercise 2
%What is the probability that if someone gets a positive test, that person is infected?

% Initialize arrays to store results
priors = 0:0.1:1;  % Prior infection rates from 0 to 1 in steps of 0.1
posterior_probs = zeros(size(priors));

% Calculate P(B) using the false positive rate
fp_rate = 0.05;
p_b = @(p_a) 1 * p_a + fp_rate * (1 - p_a);

% Calculate P(A | B) for each prior
for i = 1:length(priors)
    p_a = priors(i);
    p_b_given_a = 1;  % Probability of positive test given infection (1 because FN = 0) 
    p_a_given_b = p_b_given_a * p_a / p_b(p_a);
    posterior_probs(i) = p_a_given_b;
end

% Display results
disp('Prior (P(A))    Posterior (P(A | B))');
disp([priors; posterior_probs]);
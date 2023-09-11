% Kristen Park 
% September 10, 2023
% NGG 6050 
% Homework #3

clear
close all 

%% Exercise 1
% Assume that there are 10 quanta available in a nerve terminal, and for a 
% given release event each is released with a probability of 0.2. For one 
% such event, what is the probability that 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, or 
% 10 quanta will be released?

fprintf('\nExercise 1 \n');
n = 10;
p = 0.2; 
outcomes = 0:10; 
probabilities = binopdf(outcomes, n, p);

for i = 1:length(outcomes)
    fprintf('P(X = %d) = %.4f\n', outcomes(i), probabilities(i));
end

%% Exercise 2
% Let's say you know that a given nerve terminal contains exactly 14 quanta 
% available for release. You have read in the literature that the release 
% probability of these quanta is low, say 0.1. To assess whether this value 
% is reasonable, you run a simple experiment: activate the nerve and 
% measure the number of quanta that are released. The result is 8 quanta. 
% What is the probability that you would get this result (8 quanta) if the 
% true probability of release really was 0.1? What about if the true 
% release probability was much higher; say, 0.7? What about for each decile 
% of release probability (0.1, 0.2, ... 1.0)? Which value of release 
% probability did you determine to be the most probable, given your 
% measurement?

fprintf('\nExercise 2\n');
n = 14;
outcome = 8;
p = 0.1; 
prob = binopdf(outcome, n, p); 
fprintf('If the true probability of release was %.1f, the probability of getting this result would be %f\n\n', p, prob);

p = 0.7;
prob = binopdf(outcome, n, p);
fprintf('If the true probability of release was %.1f, the probability of getting this result would be %f\n\n', p, prob);

p = 0.1:0.1:1;
prob = binopdf(outcome,n, p);
for i = 1:length(p)
    fprintf('If the true probability of release was %.1f, the probability of getting this result would be %f\n', p(i), prob(i));
end

maxProb = max(prob);
mostProbableInd = find(prob == maxProb);
fprintf('\nThe most probable release probability is %.1f\n\n', p(mostProbableInd));

%% Exercise 3
% Not feeling convinced by your single experiment (good scientist!), you 
% repeat it under identical conditions. This time you measure 5 quanta that 
% were released. Your sample size has now doubled, to two measurements. 
% You now want to take into account both measurements when you assess the 
% likelihoods of different possible values of the underlying release 
% probability. To do so, assume that the two measurements in this sample 
% are independent of one another; that is, the value of each result had no 
% bearing on the other. In this case, the total likelihood is simply the 
% product of the likelihoods associated with each separate measurement. 
% It is also typical to compute the logarithm of each likelihood and 
% take their sum, which is often more convenient. What are the values of 
% the total likelihood and total log-likelihood in this example, if we 
% assume that the true release probability is 0.1?

fprintf('\nExercise 3\n');
n = 14;
outcome1 = 8;
outcome2 = 5;
p = 0.1; 
prob1 = binopdf(outcome1, n, p); 
prob2 = binopdf(outcome2, n, p);
totalLikelihood = prob1 * prob2; 
logLikelihood = log(prob1) + log(prob2);
fprintf('The total likelihood is %d\n', totalLikelihood);
fprintf('The log-likelihood is %f\n\n', logLikelihood);

% Of course, knowing those values of the likelihood and log-likelihood is 
% not particularly useful until you can compare them to the values computed
% for other possible values for the release probability, so you can determine 
% which value of release probability is most likely, given the data. 
% Therefore, compute the full likelihood and log-likelihood functions using 
% deciles of release probability between 0 and 1. What is the maximum value? 
% Can you improve your estimate by computing the functions at a higher 
% resolution? How does the estimate improve as you increase the sample 
% size?

p = 0.1:0.1:1;
prob1 = binopdf(outcome1, n, p);
prob2 = binopdf(outcome2, n, p); 
totalLikelihood = prob1 .* prob2; 
logLikelihood = log(prob1) + log(prob2);

for i = 1:length(p)
    fprintf('If the true probability of release was %.1f, the total likelihood would be %f and the log-likelihood would be %f\n', p(i), totalLikelihood(i), logLikelihood(i));
end

maxLikelihood = max(totalLikelihood);
fprintf('\nThe max total likelihood is %f\n\n', maxLikelihood);

p = 0.1:0.01:1;
prob1 = binopdf(outcome1, n, p);
prob2 = binopdf(outcome2, n, p); 
totalLikelihood = prob1 .* prob2; 
maxLikelihood = max(totalLikelihood);
fprintf('When using steps of 0.01, the max total likelihood is %f\n\n', maxLikelihood);

p = 0.1:0.001:1;
prob1 = binopdf(outcome1, n, p);
prob2 = binopdf(outcome2, n, p); 
totalLikelihood = prob1 .* prob2; 
maxLikelihood = max(totalLikelihood);
fprintf('When using steps of 0.001, the max total likelihood is %f\n\n', maxLikelihood);

p = 0.1:0.0001:1;
prob1 = binopdf(outcome1, n, p);
prob2 = binopdf(outcome2, n, p); 
totalLikelihood = prob1 .* prob2; 
maxLikelihood = max(totalLikelihood);
bestp = p(find(totalLikelihood == maxLikelihood)); %0.4643 most probable release probability
fprintf('When using steps of 0.0001, the max total likelihood is %f at a %.4f release probability\n\n', maxLikelihood,bestp);

% The estimate increases when computing the functions at a
% higher resolution (maxes out around 0.024056), the estimate would only
% further increase with a larger sample size. 


%% Exercise 4
% What is the most likely value of p (which we typically refer to as p^,
% which is pronounced as "p-hat" and represents the maximum-likelihood 
% estimate of a parameter in the population given our sample with a 
% resolution of 0.01?

fprintf('\nExercise 4\n');
n = 14;
observed_counts = [0, 0, 3, 10, 19, 26, 16, 16, 5, 5, 0, 0, 0, 0, 0];
release_probabilities = 0.1:0.01:1;
likelihoods = zeros(size(release_probabilities));

for i = 1:length(release_probabilities)
    p = release_probabilities(i); 
    likelihood = prod(binopdf(0:14, n, p) .^ observed_counts);
    likelihoods(i) = likelihood;
end

[max_likelihood, max_likelihood_index] = max(likelihoods);
most_likely_p = release_probabilities(max_likelihood_index);
disp(['Most likely value of p (p-hat): ' num2str(most_likely_p)]);

%% Exercise 5
% Let's say that you have run an exhaustive set of experiments on this 
% synapse and have determined that the true release probability is 0.3 
% (within some very small tolerance). Now you want to test whether changing
% the temperature of the preparation affects the release probability. 
% So you change the temperature, perform the experiment, and measure 7 
% quantal events for the same 14 available quanta. Compute p-hat. Standard 
% statistical inference now asks the question, what is the probability that 
% you would have obtained that measurement given a Null Hypothesis of no 
% effect? In this case, no effect corresponds to an unchanged value of the 
% true release probability (i.e., its value remained at 0.3 even with the 
% temperature change). What is the probability that you would have gotten 
% that measurement if your Null Hypothesis were true? Can you conclude that 
% temperature had an effect?

fprintf('\nExercise 5\n');

n = 14;               
observed = 7;   

observed_counts = zeros(1,15); %from 0 to 14
observed_counts(8) = 1; %actually for 7 because index starts at 0
release_probabilities = 0.1:0.01:1;
likelihoods = zeros(size(release_probabilities));

for i = 1:length(release_probabilities)
    p = release_probabilities(i); 
    likelihood = prod(binopdf(0:14, n, p) .^ observed_counts);
    likelihoods(i) = likelihood;
end

[max_likelihood, max_likelihood_index] = max(likelihoods);
most_likely_p = release_probabilities(max_likelihood_index);
disp(['Most likely value of p (p-hat): ' num2str(most_likely_p)]);

%null hypothesis = no effect with change in temperature 
p = 0.3;
null_prob = binopdf(observed, n, p);
disp(['The probability of obtaining this measurement if the null hypothesis were true is: ', num2str(null_prob)])

% The p-value is greater than 0.05, suggesting that one cannot reject the
% null hypothesis. The new temperature does not have a significant effect
% on release probability. 



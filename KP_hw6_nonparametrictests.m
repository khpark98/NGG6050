% Kristen Park
% 9/19/23
% Homework #6
% Simple Non-Parametric Tests


% Scenario 1: You are a behavioral geologist measuring the reaction time of 
% rocks in response to a tone. Specifically, you want to compare the median 
% reaction time of geodes, to that of limestone. You recruit 20 rocks per 
% group, and run your reaction-time experiment. What test would you use to 
% compare median reaction times between geodes and limestone, and why?

% You would use the Mann-Whitney test to compare the median reaction times of two independent groups. 


% Scenario 2: You are a brilliant scientist working at a biotech firm 
% developing a vaccine that reverses aging. Wow! To test the efficacy of 
% the vaccine, you recruit 50 people, give them a course of your vaccine, 
% and measure their age with a (very) special scale before and after 
% treatment. You want to start by refuting the simple that that the 
% participants' measured ages are not changed by the treatment. What test 
% do you use and why?

% You would use the Wilcoxon Signed-Rank test because you are dealing with non-normal, paired data. 


% Scenario 3: You are a neuroeconomist and believe you have developed a 
% wearable device that can change consumer preferences about a given 
% product. To test your device, you present product X to a group of 40 
% individuals, and ask them to fill out a survery assessing how much they 
% like the product (larger score means they like it more). Then, you have 
% the individuals wear the device, present product X, and assess how much 
% they like of the product. You want to know if the device reliably 
% increases, decreases, or does not affect their liking of product X. What 
% test would you use and why? What result would indicate that their liking 
% has increased?

% You would use the Sign test to assess for systematic direction of 
% treatment effect. A positive test statistic with a significant p-value 
% would indicate that their liking has increased. 
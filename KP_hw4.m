% Kristen Park 
% September 13, 2023
% NGG 6050 
% Homework #3

clear
close all 

% Compute confidence/credible intervals based on the four methods above for
% simulated data sampled from a population that is Gaussian distributed 
% with mean =10 and standard deviation =2, for n=5, 10, 20, 40, 80, 160, 
% 1000 at a 95% confidence level.

mu = 10; %mean
sigma = 2; %stdev 
n = [5, 10, 20, 40, 80, 160, 1000]; %n 
alpha = 0.95; 


for i = 1:length(n)
    % Generate simulated data
    data = normrnd(mu, sigma, [1, n(i)]);
    sampleMean = mean(data);
    fprintf('\nn = %i, mean = %.1f\n', n(i), sampleMean);

    %% 1. The simple, analytic approach with large n and/or known standard deviation. 
    %1a with given population stdev (sigma) 
    sem = sigma / sqrt(n(i));
    zscore = norminv(1-alpha/2); 

    lower = sampleMean - (zscore * sem);
    upper = sampleMean + (zscore * sem);
    
    fprintf('1a (with sigma) : CI=[%.3f, %.3f]\n', lower, upper);

    %1b with sample stdev (s) 
    sem = std(data) / sqrt(n(i));
    zscore = norminv(1-alpha/2); 

    lower = sampleMean - (zscore * sem);
    upper = sampleMean + (zscore * sem);

    fprintf('1b (with sample s) : CI=[%.3f, %.3f]\n', lower, upper);

    %% 2. The simple, analytic approach with small n and unknown population standard deviation
    sem = std(data) / sqrt(n(i)); %use sample stdev when using t score
    tscore = tinv(1-alpha/2, n(i)-1); %n-1 degrees of freedome

    lower = sampleMean - (tscore * sem);
    upper = sampleMean + (tscore * sem);

    fprintf('2 : CI=[%.3f, %.3f]\n', lower, upper);    

    %% 3. Bootstrapped confidence intervals
    numBootstraps = 1000; 
    muStar = bootstrp(numBootstraps, @mean, data); %means of bootstrapped population

    lower = prctile(muStar, 100 * (1 - alpha) / 2);
    upper = prctile(muStar, 100 * (alpha + (1 - alpha) / 2));

    fprintf('3 : CI=[%.3f, %.3f]\n', lower, upper);

    %% 4. Bayesian credible intervals 


end

 
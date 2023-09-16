% Kristen Park 
% September 8, 2023
% NGG 6050 
% Homework #2

% "Basal forebrain circuit for sleep-wake control"
% Xu et al. 2015
% https://www.nature.com/articles/nn.4143

%% Gaussian Simulation of Spike Latency 

% Parameters
mu = 2.1;     % Mean
sigma = 1.1;  % Standard Deviation
n = 85; %85 neurons 

% Generate random numbers from a standard normal distribution
rand_values = randn(1, num_samples);
samples = mu + sigma * rand_values;

% Plot histogram for a given number of bins (using trapz to approximate pdf)
nbins = 20;
[counts, edges] = histcounts(samples, nbins);
xaxis = (edges(1:end-1) + edges(2:end))/2;
n_pdf = counts / trapz(xaxis, counts);
bar(xaxis, n_pdf);

% Show theoretical pdf in red
hold on;
plot(xaxis, normpdf(xaxis, mu, sigma), 'r-', 'LineWidth', 2);
hold off;

% Labels, etc.
title(sprintf('Gaussian pdf, mu=%.2f, sigma=%.2f', mu, sigma));
xlabel('Value');
ylabel('Probability');
legend('Simulated', 'Theoretical');
grid on;
% Parameters
lambda = 10; % Average number of base stations (Poisson parameter)
num_sim = 10000; % Number of SINR points to simulate

% Generating Poisson-distributed number of base stations
numBS = poissrnd(lambda, 1, num_sim);

% Parameters for SINR calculation
signalPower = 1; % Signal power at receiver
interferencePower = 0.5; % Average interference power at receiver
noisePower = 0.1; % Noise power at receiver

% Calculating SINR for each point
SINR = signalPower ./ (interferencePower * (numBS-1) + noisePower);

% Sorting SINR values for CDF
sortedSINR = sort(SINR);

% Calculating CDF
cdf = (1:num_sim) / num_sim;

% Plotting CDF
figure;
plot(sortedSINR, cdf, 'b-', 'LineWidth', 2);
grid on;
xlabel('SINR');
ylabel('CDF');
title('CDF of SINR with Poisson-distributed Base Stations');

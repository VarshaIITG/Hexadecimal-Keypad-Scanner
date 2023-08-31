clear all;
close all;

% Parameters
lambda = 0.01; % Base station density (Poisson parameter)
numsim = 100000; % Number of simulations for each SINR threshold
threshold_SINR_dB = -10:2:20; % Threshold SINR values in dB
threshold_SINR_lin = 10.^(threshold_SINR_dB / 10); % Convert to linear scale

% Initialize outage probability vector
outage_prob = zeros(size(threshold_SINR_dB));

% Simulation loop for different threshold SINR values
for i = 1:length(threshold_SINR_dB)
    out_num = 0;

    for N = 1:numsim
        % Generate Poisson distributed base station locations
        numBaseStations = poissrnd(lambda);
        bs_locations = rand(numBaseStations, 2); % Uniform distribution [0, 1]

        % Calculating distances and received powers
        dis = sqrt(sum(bs_locations.^2, 2));
        rx_powers = 1 ./ dis.^2;

        % Finding the nearest base station to the user
        [~, nearest_bs_indx] = min(dis);
        
        % Calculate interference from other base stations
        interference = sum(rx_powers) - rx_powers(nearest_bs_indx);

        % Calculate SINR
        sinr = rx_powers(nearest_bs_indx) / (interference + 1);

        % Check if SINR meets the threshold
        if sinr < threshold_SINR_lin(i)
            out_num = out_num + 1;
        end
    end

    % Calculate outage probability
    outage_prob(i) = out_num / numsim;
end

% Plot outage probability vs. threshold SINR
figure;
semilogy(threshold_SINR_dB, outage_prob, 'bo-');
grid on;
xlabel('Threshold SINR (dB)');
ylabel('Outage Probability');
title('Outage Probability vs. Threshold SINR in Poisson cellular Network');

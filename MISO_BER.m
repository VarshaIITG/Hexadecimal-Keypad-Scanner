clear all;
close all;

% Parameters
N = 100000; % Number of bits to simulate
snr_dB = 0:1:20; % Range of SNR values in dB
snr_lin = 10.^(snr_dB/10); % Convert SNR values to linear scale
Tx = [1,2,4,8]; % Different values of N (number of transmit antennas)

% Generate random bits
bits = randi([0, 1], 1, N);

% Modulation(BPSK)
modulatedBits = 2 * bits - 1;

% Initializing BER matrix
ber_matrix = zeros(length(Tx), length(snr_dB));

% Simulation loop for different values of NT
for nTx = 1:length(Tx)
    NT = Tx(nTx);
    
    % Generate Rayleigh fading channel coefficients
    H = (randn(1, NT) + 1i * randn(1, NT)) / sqrt(2);
    
    % BER vector for this N
    ber = zeros(size(snr_dB));
    
    % Simulation loop for different SNR values
    for i = 1:length(snr_dB)
        % Add complex Gaussian noise to the transmitted signal
        noise = sqrt(1 / (2 * snr_lin(i))) * (randn(1, N) + 1i * randn(1, N));
        rx_signal = (H'*diag(sqrt(snr_lin(i)))).*( modulatedBits) + noise;



        % Equalization (using matched filter)
        demod_Bits = sum(conj(H') .* rx_signal, 1);

        % Hard decision decoding
        decodedBits = real(demod_Bits) > 0;

        % Calculating BER
        ber(i) = sum(decodedBits ~= bits) / N;
    end
    
    % Storing BER values for this NT
    ber_matrix(nTx, :) = ber;
end

% Plotting BER vs. SNR for different values of N
figure;
semilogy(snr_dB, ber_matrix(1, :), 'bo-', 'DisplayName', 'N=1');
hold on;
semilogy(snr_dB, ber_matrix(2, :), 'go-', 'DisplayName', 'N=2');
semilogy(snr_dB, ber_matrix(3, :), 'ro-', 'DisplayName', 'N=4');
semilogy(snr_dB, ber_matrix(4, :), 'co-', 'DisplayName', 'N=8');
grid on;
hold off;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR for N x 1 MISO Rayleigh Fading Channel');
legend('Location', 'best');

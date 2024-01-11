clc;
clear all;
close all;

% Parameters

N = 2;          % Number of antennas per Access Point
rho = 1;        % Transmit power (W)
tau_p = 10;     % Length of pilot sequence
K_values = 10:10:100; % Range of values for K (number of USERS)
sigma_degree = 10; % Sigma in degrees
sigma_radian = deg2rad(sigma_degree);  % Convert sigma from degrees to radians
d = 20; % Distance parameter d (you can adjust this as needed)
M = 10;

% Additional Parameters
nf_dB = 9;      % Noise figure (dB)
B = 20e6;       % System bandwidth (Hz)
eta_k = 1;      % Power control coefficient
P_tc_m = 0.2;   % Transmit power of tc (W)
P_0_m = 0.825;  % Constant power P_0_m (W)
P_bt_m = 0.25 / 1e9; % Transmit power per Gbits per second (W)
kB = 1.381e-23; % Boltzmann constant (J/K)
T0 = 290;       % Reference temperature (K)

% Additional Parameters for Monte Carlo Simulation
numsim = 10; % Number of Monte Carlo simulations

% Simulation
EE_values = zeros(numsim, length(K_values));
EE_value = zeros(numsim, length(K_values));

for sim = 1:numsim
    % Initialize EE_sum for this simulation
    EE_sum = 0;
    SE_k_sum = 0;

    for K_idx = 1:length(K_values)
        K = K_values(K_idx);

        % Initialize total power
        P_total = 0;

        for sim_num = 1:numsim
           for k = 1:K
           
            %Generating channel 
            h=1/sqrt(2).*[randn(N,1)+1j*randn(N,1)];
            Beta=randn(N,1);
            g=sqrt(Beta).*h;

            % Uplink Training
            tau_up=M;
            ps=randi([0,1],M,1);
            w_p=1/sqrt(2).*[randn(N,1)+1j*randn(N,1)];
            Rho=20;

            %  channel estimation
            w_ud=1/sqrt(2).*[randn(N,1)+1j*randn(N,1)];
            g_Hat=g+w_ud/sqrt(tau_up.*Rho); 
            % Calculate large scale fading coefficient betamk using the model given in the paper
            d_mk = sqrt(rand^2 + rand^2) * 1e3; % Convert to meters
            F_mk = sqrt(100) * randn; % Mean 0, variance 100 (dB)
            beta_mk_dB = -34.53 - 38 * log10(d_mk / 1) + F_mk;
            beta_mk = 10^(beta_mk_dB / 20); % Convert to linear scale

            % Generate spatial correlation matrix Rmk
            phi = -pi + (2*pi) * rand(1, N);
            corr_matrix = zeros(N, N);

            for l = 1:N
                for n = 1:N
                    A = 2 * pi * d * (l - n) * cos(phi(n)) * sigma_radian;
                    term1 = exp(2 * pi * 1i * d * (l - n) * sin(phi(n)));
                    term2 = exp((-sigma_radian^2 / 2) * (2 * pi * d * (l - n) * cos(phi(n)))^2);
                    term3 = qfunc((-20 * sigma_radian - A) / sigma_radian) - qfunc((20 * sigma_radian - A) / sigma_radian);
                    corr_matrix(l, n) = beta_mk * (term1 * term2 * term3);
                end
            end

            % Ensure the matrix is Hermitian (conjugate symmetric)
            Rmk = corr_matrix + conj(corr_matrix') - diag(diag(corr_matrix));

            % Calculate SE for each K
            DSk = sqrt(rho) * mean(diag(Rmk));
            BUk = sqrt(rho) * (trace(Rmk^2) - trace(Rmk)^2 / N);
            UIkk = sqrt(rho) * trace(Rmk);
            NIk = mean(diag(Rmk));

            SE_k = (abs(DSk)^2) / (abs(BUk)^2 + abs(UIkk)^2 + abs(NIk)^2);
            ss = log2(1 + SE_k);
            
            % Calculate FH power
            P_fh_m = P_0_m + B * M * ss * P_bt_m;

            % Calculate total power
            P_total = P_total + rho * M + N * K * P_tc_m + K * P_fh_m;
            
            % Calculate EE for each user
            EE_k = B * M * ss / P_total;

            % Accumulate EE for each user
            EE_sum = EE_sum + EE_k;
        end
        

            % Accumulate SE for each user
            SE_k_sum = SE_k_sum + ss;
        end

        % Average EE over all users and APs and simulations
        EE_values(sim, K_idx) = EE_sum / (M * K * numsim);
        avg_ss = SE_k_sum / (M * numsim);

        % Calculate FH power for average SS
        P_fh_avg = P_0_m + B * M * avg_ss * P_bt_m;

        % Calculate total power for average SS
        P_total_avg = rho * M + N * K * P_tc_m + K * P_fh_avg;

        % Calculate EE for average SS
        EE_value(sim, K_idx) = B * M * avg_ss / P_total_avg;
    end
end

% Average results over Monte Carlo simulations
average_EE_values = mean(EE_values, 1);
average_EE_value = mean(EE_value, 1);

% Plot Results
figure;
set(0,'DefaultAxesFontName', 'Times New Roman');
plot(K_values, average_EE_value, 's-', 'LineWidth', 2, 'DisplayName', 'Average EE');
xlabel('Number of Users');
ylabel('Energy Efficiency (bit/J)');
% title('Energy Efficiency vs. Number of Users in Cell-Free Massive MIMO');
legend('show');

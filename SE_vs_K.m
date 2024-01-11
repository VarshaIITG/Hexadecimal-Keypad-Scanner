clc;
clear all;
close all;

% Parameters
M = 10;         % Number of Access Points
N = 2;          % Number of antennas per Access Point
rho = 1;        % Normalized uplink SNR
    
tau_c = 100;    % Coherence time
K_values = 5:5:50; % Range of values for K (number of USERS)

sigma_degree = 10; % Sigma in degrees
sigma_radian = deg2rad(sigma_degree);  % Convert sigma from degrees to radians
d = 20;         % Distance parameter d (we can adjust this as needed according to distance in antenna array)

numsim = 100; % Number of Monte Carlo simulations

% Simulation

ss_values = zeros(size(K_values));

for k_idx = 1:length(K_values)
    K = K_values(k_idx);
    tau_p = K; % Length of pilot sequence
   
    % Initialize SE_sum for this simulation
    ds=0;
    bu=0;
    ui=0;
    ni=0;
    SINR_sum = 0;
    
    for sim = 1:numsim
        for m = 1:M
            %Generating channel 
            h=1/sqrt(2).*[randn(N,1)+1j*randn(N,1)];
            Beta=randn(N,1);
            g=sqrt(Beta).*h;

            % Uplink Training
            tau_up=K;
            ps=randi([0,1],K,1);
            w_p=1/sqrt(2).*[randn(N,1)+1j*randn(N,1)];
            Rho=20;

            %  channel estimation
            w_ud=1/sqrt(2).*[randn(N,1)+1j*randn(N,1)];
            g_Hat=g+w_ud/sqrt(tau_up.*Rho);            

            % Generating large scale fading coefficient betamk using the model given in the paper
            d_mk = sqrt(rand^2 + rand^2) * 1e3; % Convert to meters
            F_mk = sqrt(100) * randn; % Mean 0, variance 100 (dB)
            beta_mk_dB = -34.53 - 38 * log10(d_mk / 1) + F_mk;
            beta_mk = 10^(beta_mk_dB / 20); % Convert to linear scale

            % Generate spatial correlation matrix R_mk
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
            Rmk = corr_matrix ;

            % Calculate SE for each K
            DSk = ds+sqrt(rho) * mean(diag(Rmk));
            BUk = bu+sqrt(rho) * (trace(Rmk^2) - trace(Rmk)^2 / N);
            UIkk = ui+sqrt(rho) * trace(Rmk);
            NIk = ni +mean(diag(Rmk));

            SINR_k = (abs(DSk)^2) / (abs(BUk)^2 + abs(UIkk)^2 + abs(NIk)^2);

            % Accumulate SE for each user
            SINR_sum = SINR_sum + SINR_k;
        end
    end

    % Average SE over all users and APs
    ss_values(k_idx) = (1-tau_p/tau_c).*log2(1 + SINR_sum / (K * numsim));
end

% Plot Results
figure;

plot(K_values, ss_values, 's-', 'LineWidth', 2, 'DisplayName', 'ss');
xlabel('Number of Users (K)');
ylabel('Spectral Efficiency (SE)');
% title('Spectral Efficiency vs. Number of Users in Cell-Free Massive MIMO');
legend('show');
grid on;

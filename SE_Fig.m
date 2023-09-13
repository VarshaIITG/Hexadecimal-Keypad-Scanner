clc;
M = 10;         % Number of Accesss points
N = 4;
numsim = 1e1;
K_values = 5:5:50; % Range of values for K (number of USERS)
sigma_degree = 30; % sigma in degrees
sigma_radian = deg2rad(sigma_degree);  % Convert sigma from degrees to radians
d = 20;         % Distance parameter d

% Initialize SE_corr array
SE_corr = zeros(1, length(K_values));

% Iterate over different values of users
for ku = 1:length(K_values)
    K = K_values(ku);
    SE_k = 0;
    for k=1:K

    for sim = 1:numsim
        dS_k = 0;
        bU_k = 0;
        nI_k = 0;
        SINR_k = 0;
        

        for m = 1:M
            % Generate spatial correlation matrix R_mk
            phi = -pi + (2*pi) * rand(N, N);
            corr_matrix = zeros(N, N);
            

            for l = 1:N
                for n = 1:N
                     A=2*pi*d*(l-n)*cos(phi)*sigma_radian;
                     term1 = exp(2*pi*1j*d*(l - n)*sin(phi(l, n)));
                     term2 = exp((-sigma_radian^2/2)*(2*pi*d*(l - n)*cos(phi(l, n)))^2);
                     term3 = qfunc((-20*sigma_radian-A)/sigma_radian)-qfunc((20*sigma_radian-A)/sigma_radian);
                     corr_matrix(l, n) = (term1 * term2*term3);
                 end
             end
             
             % Ensure the matrix is Hermitian (conjugate symmetric)
             R_mk=corr_matrix;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ;
            %corr_matrix=eye(N);

            % Calculate large scale fading coefficient betamk using the model given in paper
            d_mk = sqrt(rand^2 + rand^2) * 1e3; % Convert to meters
            F_mk = sqrt(100) * randn; % Mean 0, variance 100 (dB)
            beta_mk_dB = -34.53 - 38 * log10(d_mk / 1) + F_mk;
            beta_mk = 10^(beta_mk_dB / 20); % Convert to linear scale

            % Generate a complex Gaussian vector h_mk with mean 0 and covariance matrix as the identity matrix
            h_mk = sqrt(1/2) * (randn(N, 1) + 1i * randn(N, 1));
            w_mk = sqrt(1/2) * (randn(N, 1) + 1i * randn(N, 1));

            % Compute the channel vector
            g_mk = sqrt(beta_mk) * (sqrtm(R_mk) * h_mk);
            gcap_mk = g_mk + w_mk;

            dS_k = dS_k + g_mk' * g_mk;
            bU_k = bU_k + g_mk' * g_mk + g_mk' * w_mk;
            nI_k = nI_k + gcap_mk' * w_mk;
        end

        % Calculate the expectation by dividing the sum by M
        DS_k = dS_k / M;
        BU_k = bU_k - DS_k;
        NI_k = nI_k;
        SINR_k = SINR_k + ((abs(DS_k))^2 / ((abs(BU_k))^2 + (abs(NI_k))^2));
        SE_k = SE_k + log2(1 + SINR_k);
    end
    
    end
    % Average SE over simulations and number of users K
    SE_corr(ku) = SE_k / K;
end
% Plot the results
plot(K_values, SE_corr, 'b-o');
xlabel('Number of users(k)');
ylabel('Average Spectral Efficiency (SE)');
legend('correlated Fading');
title('Average SE vs. Number of APs');
grid on;

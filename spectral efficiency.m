% Constants and parameters
K = 10;         % Number of users
M_array = 1:10:100; % Range of values for M (number of APs)
sigma_degree = 30; % sigma in degrees
% Convert sigma from degrees to radians
sigma_radian = deg2rad(sigma_degree);
d = 20;         % Distance parameter d
Fmk = normrnd(0, 10);

% Iterate over different values of M
for m_index = 1:length(M_array)
    M = M_array(m_index); % Number of APs
    a_m = zeros(M, 1); % Initialize array for 'a'
    % Generating random values for phi between -pi and pi
    phi = -pi + (2*pi) * rand(M, M); % Generate phi with dimensions (M, M)
    % Initialize the correlation matrix
    corr_matrix = zeros(M, M);
    
    for k = 1:K
        for m = 1:M
            % Calculate correlation matrix elements
            for l = 1:M
                for n = 1:M
                    term1 = exp(2*pi*1j*d*(l - n)*sin(phi(l, n)));
                    term2 = exp((-sigma_radian^2/2)*(2*pi*d*(l - n)*cos(phi(l, n)))^2);
                    corr_matrix(l, n) = (term1 * term2).*(term1 * term2);
                end
            end
            R_mk = corr_matrix;

            AP_coordinates = rand(M, 2); % M random (x, y) positions for APs
            user_coordinates = rand(K, 2); % K random (x, y) positions for users

            % Initialize an array to store the large-scale fading coefficients betamk
            beta_mk = zeros(M, K);

            % Calculate the large-scale fading coefficient betamk for each pair
            for k = 1:K
                for m = 1:M
                    % Calculate the distance dmk between the m-th AP and k-th user
                    dmk = norm(AP_coordinates(m, :) - user_coordinates(k, :));
                    % Calculate betamk using the model given in paper
                    beta_mk(m, k) = -34.53 - 38*log10(dmk/1) + Fmk;
                end
            end
            
            a_m(m) = log2(abs(a_m(m)+(beta_mk(m, k)*trace(R_mk))).^2);
        end
    end
    
    % Calculate the average 'a' over k users for this specific M
    a(m_index) = mean(a_m);
end

% Plot the results
figure;
plot(M_array, a, 'o-');
xlabel('Number of APs (M)');
ylabel('Value of a');
title('Value of a vs. Number of APs');
grid on;

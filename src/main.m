clc;
clear;
close all; 

%% Importing data from csv file
filename = fullfile('..','data', 'group_08.csv');
opts = detectImportOptions(filename, ...
    'Delimiter', ';', ...
    'DecimalSeparator', ',');
T = readtable(filename, opts);

t       = T{:,1};
q1      = T{:,2};
h1      = T{:,3};
h2      = T{:,4};
q2_meas = T{:,5};

%% Identification of model parameters
% Noise parameters
noise_only = q2_meas(1:50);         
E_hat = mean(noise_only);           % Estimate of noise mean value
sig_hat =  var(noise_only);         % Estimate of noise variance
noise_std  = std(noise_only);
fprintf('=============== Noise parameters ==============\n'); 
fprintf('Noise mean value = %.6f\n', E_hat); 
fprintf('Noise variance = %.6f\n', sig_hat);
fprintf('Noise standard deviation = %.6f\n', noise_std);


% 'K' parameter identification (LS)
idx = 51:length(t);          
sqrt_h2 = sqrt(max(h2(idx), 0));
numerator   = sum(sqrt_h2 .* q2_meas(idx));
denominator = sum(h2(idx));
k_hat = numerator / denominator;
fprintf('=============== Model parameters ==============\n'); 
fprintf('K_hat = %.6e\n', k_hat);

% Model output for estimated k
q2_model = k_hat * sqrt_h2;

% Validation of parameter K value on figures
figure;
subplot(2,1,1)
plot(t(idx), q2_meas(idx), 'b', 'LineWidth', 1.2); hold on;
plot(t(idx), q2_model, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Outlet flow q_2');
legend('Measured q_2', 'Model: k\_hat \cdot sqrt(h_2)', 'Location', 'best');
title('Measured and modeled q2 flow in time');
subplot(2,1,2)
plot(sqrt_h2, q2_meas(idx), 'bo'); hold on;
plot(sqrt_h2, q2_model, 'r.', 'MarkerSize', 12);
grid on;
xlabel('sqrt(h_2)');
ylabel('Outlet flow q_2');
legend('Measurements', 'Model', 'Location', 'best');
title('Relationship between q_2 and sqrt(h_2)');


function plot_ekf_results(t, h1, h2, q2_meas, x_hat_history, y_pred_history, innovation_history, k_hat)

% Reconstructed outlet flow from corrected state estimate
q2_est_from_state = k_hat .* sqrt(max(x_hat_history(2,:), 0));

%% Figure 1: state estimation
figure;

subplot(2,1,1)
plot(t, h1, 'b', 'LineWidth', 1.4); hold on;
plot(t, x_hat_history(1,:)', 'r--', 'LineWidth', 1.6);
grid on;
xlabel('Time');
ylabel('h_1');
legend('True h_1', 'Estimated h_1', 'Location', 'best');
title('True and estimated h_1');

subplot(2,1,2)
plot(t, h2, 'b', 'LineWidth', 1.4); hold on;
plot(t, x_hat_history(2,:)', 'r--', 'LineWidth', 1.6);
grid on;
xlabel('Time');
ylabel('h_2');
legend('True h_2', 'Estimated h_2', 'Location', 'best');
title('True and estimated h_2');

%% Figure 2: outlet flow reconstructed from corrected state
figure;

subplot(2,1,1)
plot(t, q2_meas, 'b', 'LineWidth', 1.2); hold on;
plot(t, q2_est_from_state, 'r--', 'LineWidth', 1.6);
grid on;
xlabel('Time');
ylabel('q_2');
legend('Measured q_2', 'Estimated q_2 from corrected state', 'Location', 'best');
title('Measured and estimated outlet flow');

subplot(2,1,2)
plot(q2_est_from_state, q2_meas, 'o', ...
    'MarkerEdgeColor', [0 0.45 0.74], ...
    'MarkerFaceColor', 'none');
grid on;
xlabel('Estimated q_2');
ylabel('Measured q_2');
title('Measured vs estimated outlet flow');

%% Figure 3: predicted measurement and innovation
figure;

subplot(2,1,1)
plot(t, q2_meas, 'b', 'LineWidth', 1.2); hold on;
plot(t, y_pred_history, 'm--', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('q_2');
legend('Measured q_2', 'Predicted q_2 before correction', 'Location', 'best');
title('Predicted measurement before correction');

subplot(2,1,2)
plot(t, innovation_history, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
grid on;
xlabel('Time');
ylabel('Innovation');
title('Innovation signal');

%% Figure 4: state estimation errors
figure;

subplot(2,1,1)
plot(t, h1 - x_hat_history(1,:)', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.2);
grid on;
xlabel('Time');
ylabel('Error in h_1');
title('Estimation error for h_1');

subplot(2,1,2)
plot(t, h2 - x_hat_history(2,:)', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.2);
grid on;
xlabel('Time');
ylabel('Error in h_2');
title('Estimation error for h_2');

end
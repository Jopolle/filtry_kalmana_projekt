function plot_identification_results(t, idx, idx2, h1, h2, q2_meas, sqrt_h2, q2_model, dh2, dh2_model, dh1, dh1_model, Psi)

%% Validation of parameter K
figure;
subplot(2,1,1)
plot(t(idx), q2_meas(idx), 'b', 'LineWidth', 1.2); hold on;
plot(t(idx), q2_model, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('Outlet flow q_2');
legend('Measured q_2', 'Model: k_{hat} \cdot sqrt(h_2)', 'Location', 'best');
title('Measured and modeled q_2 flow in time');

subplot(2,1,2)
plot(sqrt_h2, q2_meas(idx), 'bo'); hold on;
plot(sqrt_h2, q2_model, 'r.', 'MarkerSize', 12);
grid on;
xlabel('sqrt(h_2)');
ylabel('Outlet flow q_2');
legend('Measurements', 'Model', 'Location', 'best');
title('Relationship between q_2 and sqrt(h_2)');

%% Validation of alpha and S2
figure;
subplot(2,1,1)
plot(t(idx2), dh2, 'b', 'LineWidth', 1.2); hold on;
plot(t(idx2), dh2_model, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('dh_2/dt');
legend('Estimated derivative dh_2/dt', 'Model: \Phi \theta_{hat}', 'Location', 'best');
title('Measured and modeled dh_2/dt in time');

subplot(2,1,2)
plot(h1(idx2)-h2(idx2), dh2, 'bo'); hold on;
plot(h1(idx2)-h2(idx2), dh2_model, 'r.', 'MarkerSize', 12);
grid on;
xlabel('h_1 - h_2');
ylabel('dh_2/dt');
legend('Estimated derivative', 'Model', 'Location', 'best');
title('Relationship between dh_2/dt and h_1-h_2');

%% Validation of S1
figure;
subplot(2,1,1)
plot(t(idx2), dh1, 'b', 'LineWidth', 1.2); hold on;
plot(t(idx2), dh1_model, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Time');
ylabel('dh_1/dt');
legend('Estimated derivative dh_1/dt', 'Model: a3_{hat} \cdot \Psi', 'Location', 'best');
title('Measured and modeled dh_1/dt in time');

subplot(2,1,2)
plot(Psi, dh1, 'bo'); hold on;
plot(Psi, dh1_model, 'r.', 'MarkerSize', 12);
grid on;
xlabel('\Psi = q_1 - \alpha_{hat}(h_1-h_2)');
ylabel('dh_1/dt');
legend('Estimated derivative', 'Model', 'Location', 'best');
title('Relationship between dh_1/dt and \Psi');
end
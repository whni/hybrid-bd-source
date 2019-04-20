%%%%% Plotting
SNR = -40:5:0;

lw = 1.5;
ms = 5;

load('Alpha-Nt32-K2-Nr16-Ns1-Lr1-Lt2-Ncls1-Nray1.mat')

rate_bd_1norm = rateMat(:, 1);
rate_limited = rateMat(:, 2);
rate_alpha = rateMat(:, 3);

figure
hold on
plot(SNR, abs(rate_alpha), 'k-v', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_1norm), 'b--', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_limited), 'go-.', 'LineWidth', lw, 'MarkerSize', ms)

hold off

legend('Hybrid BD (Analytical)', 'Hybrid BD (Simulation)', 'Limited feedback hybrid precoding [18]')
xlabel('SNR (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')
grid


load('ALpha-Nt256-K2-Nr128-Ns1-Lr1-Lt2-Ncls1-Nray1.mat')
rate_bd_1norm = rateMat(:, 1);
rate_limited = rateMat(:, 2);
rate_alpha = rateMat(:, 3);

%figure
hold on
plot(SNR, abs(rate_alpha), 'k-v', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_1norm), 'b--', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_limited), 'go-.', 'LineWidth', lw, 'MarkerSize', ms)

hold off

legend('Hybrid BD (Analytical)', 'Hybrid BD (Simulation)', 'Limited feedback hybrid precoding [18]')
xlabel('SNR (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')

grid on

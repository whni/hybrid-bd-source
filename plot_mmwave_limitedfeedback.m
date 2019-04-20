%%%%% Plotting
SNR = -40:5:0;

lw = 1.5;
ms = 5;

load('LF-ULA-Nt256-K8-Nr16-Ns1-Lr1-Lt8-Ncls8-Nray10.mat')
rate_bd = rateMat(:, 1);
rate_bd_spa1 = rateMat(:, 2);
rate_bd_1norm = rateMat(:, 3);
rate_limited = rateMat(:, 4);

figure
hold on
plot(SNR, abs(rate_bd), 'k-*', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_1norm), 'b--', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_spa1), 'r^-.', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_limited), 'go-.', 'LineWidth', lw, 'MarkerSize', ms)

hold off

legend('Traditional BD', 'Hybrid BD (1-norm)', ...
    'Sparse precoding&combining [15]', 'Limited feedback hybrid precoding [18]')
xlabel('SNR (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')
grid


load('LF-UPA-Nt256-K8-Nr16-Ns1-Lr1-Lt8-Ncls8-Nray10.mat')
rate_bd = rateMat(:, 1);
rate_bd_spa1 = rateMat(:, 2);
rate_bd_1norm = rateMat(:, 3);
rate_limited = rateMat(:, 4);

figure
hold on
plot(SNR, abs(rate_bd), 'k-*', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_1norm), 'b--', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_spa1), 'r^-.', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_limited), 'go-.', 'LineWidth', lw, 'MarkerSize', ms)

hold off

legend('Traditional BD', 'Hybrid BD (1-norm)', ...
    'Sparse precoding&combining [15]', 'Limited feedback hybrid precoding [18]')
xlabel('SNR (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')

grid

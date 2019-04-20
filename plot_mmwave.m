%%%%% Plotting
SNR = -40:5:0;
lw = 1.5;
ms = 5;

load('t256r16k8s2c2_mmwave.mat')
rate_bd_fnorm = rateMat(:, 3);
rate_bd = rateMat(:, 1);
rate_bd_spa2 = rateMat(:, 2);

load('t256r16k8s2c4_mmwave.mat')
rate_bd_spa4 = rateMat(:, 2);

figure
hold on
plot(SNR, abs(rate_bd), 'k-*', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_fnorm), 'b--', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_spa2), 'r^-.', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_spa4), 'rv-.', 'LineWidth', lw, 'MarkerSize', ms)

hold off
grid

legend('Traditional BD', 'Hybrid BD (1-norm)', ...
    'Sparse precoding&combining (2 RF chains) [15]', 'Sparse precoding&combining (4 RF chains) [15]')
xlabel('SNR (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')


load('upa_t256r16k8s2c2_mmwave.mat')
rate_bd_fnorm = rateMat(:, 3);
rate_bd = rateMat(:, 1);
rate_bd_spa2 = rateMat(:, 2);

load('upa_t256r16k8s2c4_mmwave.mat')
rate_bd_spa4 = rateMat(:, 2);

figure
hold on
plot(SNR, abs(rate_bd), 'k-*', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_fnorm), 'b--', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_spa2), 'r^-.', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_spa4), 'rv-.', 'LineWidth', lw, 'MarkerSize', ms)

hold off

legend('Traditional BD', 'Hybrid BD (1-norm)', ...
    'Sparse precoding&combining (2 RF chains) [15]', 'Sparse precoding&combining (4 RF chains) [15]')
xlabel('SNR (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')

grid

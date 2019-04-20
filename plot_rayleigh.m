%%%%% Plotting
SNR = -40:5:0;
figure
lw = 1.5;
ms = 5;

load('t256r16k8s2c2_rayleigh.mat')
rate_bd_fnorm = RAY_RATE_SET(:, 2);
rate_bd = RAY_RATE_SET(:, 1);
load('t256r16k8s2c2_ray_1norm.mat')
rate_bd_1norm = RAY_RATE_SET(:, 2);

hold on
plot(SNR, abs(rate_bd), 'k-*', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_1norm), 'b--', 'LineWidth', lw, 'MarkerSize', ms)
%plot(SNR, abs(rate_bd_fnorm), 'ro', 'LineWidth', lw, 'MarkerSize', ms)
hold off

legend('Traditional BD', 'Hybrid BD (1-norm)')
xlabel('SNR (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')



load('t64r4k8s2c2_rayleigh.mat')
rate_bd_fnorm = RAY_RATE_SET(:, 2);
rate_bd = RAY_RATE_SET(:, 1);
load('t64r4k8s2c2_ray_1norm.mat')
rate_bd_1norm = RAY_RATE_SET(:, 2);

hold on
plot(SNR, abs(rate_bd), 'k-*', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rate_bd_1norm), 'b--', 'LineWidth', lw, 'MarkerSize', ms)
%plot(SNR, abs(rate_bd_fnorm), 'ro', 'LineWidth', lw, 'MarkerSize', ms)
hold off

legend('Traditional BD', 'Hybrid BD (1-norm)')
xlabel('SNR (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')

grid

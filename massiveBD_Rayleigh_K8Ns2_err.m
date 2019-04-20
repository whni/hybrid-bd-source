% Study capacity performance of BD in massive MIMO
% Hybrid precoding VS full-complexity BD
% By Weiheng Ni, UVic, Oct. 27, 2013
% Modified July 16, 2014

tic; clear all; clc;

% =============================================
% ============ Parameter settings =============
% =============================================
Nt = 256;
Nr = 16;

SNR = 0;
K = 8;      % User number
Ns = 4;     % #streams per user
Lr = 4;     % #chains per user
Lt = Lr*K;  % #chain at BS

% randn('state', 3);
channNum = 50;  % For each user
display('Channel Loading ...')
genH = zeros(Nr, Nt, channNum*K);
for ch = 1:channNum*K
    genH(:, :, ch) = (randn(Nr, Nt) + j*randn(Nr, Nt)) / sqrt(2);
end
display('Successful!');

SNR_err = -40: 2 : -6;
nSNR = length(SNR_err);

j = sqrt(-1);

rateBD = zeros(nSNR, 1);        % full-complexity BD
ratephaseBD_dft = zeros(nSNR, 1);   % Phase + BB-BD
rateBD_err = zeros(nSNR, 1);        % full-complexity BD with err
ratephaseBD_dft_err = zeros(nSNR, 1);   % Phase + BB-BD with err

for isnr = 1 : nSNR
    P = 10^(SNR/10);
    err = 10^(SNR_err(isnr)/10);
    for ichannel = 1 : channNum   
        
        % ============= channel generation =====================
        H = zeros(K*Nr, Nt);
        Hcell = cell(K, 1);
        Hest = zeros(K*Nr, Nt);
        Herr =  err * (randn(K*Nr, Nt) + j*randn(K*Nr, Nt)) / sqrt(2);
        for k = 1:K
            Hcell{k} = genH(:, :, ichannel + channNum * (k-1));
            H((k-1)*Nr + (1:Nr), 1:Nt) = Hcell{k};
            Hest((k-1)*Nr + (1:Nr), 1:Nt) = Hcell{k};
            Hcell{k} = genH(:, :, ichannel + channNum * (k-1));
            H((k-1)*Nr + (1:Nr), 1:Nt) = Hcell{k};
        end
        
        Hest = Hest + Herr;     % add estimation error
        Rs = P/(K*Ns)*eye(K*Ns);
      
        % ================= BD precoding =====================
        [TBD, WBD] = CalPrecoderBD(H, K, Ns);% BD precoding & combining
        HH = H*TBD;
        WW = zeros(K*Ns, K*Nr);
        for k = 1:K
            WW((k-1)*Ns + (1:Ns), (k-1)*Nr + (1:Nr)) = WBD((k-1)*Ns + (1:Ns), 1:Nr);
        end
        
        vec_power_gain = diag(abs(WW*HH)).^2; 
        vec_base_level = ones(K*Ns, 1) ./ vec_power_gain;
        vec_power_alloc = water_filling(P, vec_base_level);
        Fwf = diag(sqrt(K*Ns*vec_power_alloc/P));
        TBDwf = TBD*Fwf;
        rateBD(isnr) = rateBD(isnr) + calRateMU(H, TBDwf, WBD, Rs, Ns*ones(K,1));

        
        % ================= BD precoding =====================
        [TBDest, WBDest] = CalPrecoderBD(Hest, K, Ns);% BD precoding & combining
        rateBD_err(isnr) = rateBD_err(isnr) + calRateMU(H, TBDest, WBDest, Rs, Ns*ones(K,1));
        
        
        % ==================phase + BD (DFT)====================
        Wff = zeros(Nr, K*Lr);
        
        Heq = zeros(K*Lr, Nt);
        for k = 1:K
            Wff(:, (k-1)*Lr + (1:Lr)) = select_wf_bases(Hcell{k}, 1, 1/sqrt(Nr) * dftmtx(Nr), Lr);
            Heq((k-1)*Lr + (1:Lr), :) = Wff(:, (k-1)*Lr + (1:Lr))' * Hcell{k};
        end

        Trf = exp(-j*angle(Heq)') / sqrt(Nt);
        Hbb = Heq * Trf;
        [TBD, WBD] = CalPrecoderBD(Hbb, K, Ns);
        
        HH = Hbb*TBD;
        WW = zeros(K*Ns, K*Lr);
        for k = 1:K
            WW((k-1)*Ns + (1:Ns), (k-1)*Lr + (1:Lr)) = WBD((k-1)*Ns + (1:Ns), :);
        end

        vec_power_gain = diag(abs(WW*HH)).^2; 
        vec_base_level = ones(K*Ns, 1) ./ vec_power_gain;
        vec_power_alloc = water_filling(P, vec_base_level);
        Fwf = diag(sqrt(K*Ns*vec_power_alloc/P));
        
        WphaseBD = zeros(K*Ns, Nr);
        for k = 1:K
            WphaseBD((k-1)*Ns + (1:Ns), :) = WBD((k-1)*Ns + (1:Ns), :) * (Wff(:, (k-1)*Lr + (1:Lr))');
        end

        TBD = TBD*sqrt(K*Ns)/norm(Trf*TBD*Fwf, 'fro');
        TphaseBD = Trf*TBD*Fwf;
        
        ratephaseBD_dft(isnr) = ratephaseBD_dft(isnr) + calRateMU(H, TphaseBD, WphaseBD, Rs, Ns*ones(K,1));    
        
        % =================== phase BD with error ==================
        Wff_est = zeros(Nr, K*Lr);
        Heq_est = zeros(K*Lr, Nt);
        
        for k = 1:K
            Wff_est(:, (k-1)*Lr + (1:Lr)) = select_wf_bases(Hest((k-1)*Nr + (1:Nr), 1:Nt), 1, 1/sqrt(Nr) * dftmtx(Nr), Lr);
            Heq_est((k-1)*Lr + (1:Lr), :) = Wff_est(:, (k-1)*Lr + (1:Lr))' * Hest((k-1)*Nr + (1:Nr), 1:Nt);
        end

        Trf_est = exp(-j*angle(Heq_est)') / sqrt(Nt);
        Hbb_est = Heq_est * Trf_est;
        [TBD_est, WBD_est] = CalPrecoderBD(Hbb_est, K, Ns);
    
        WphaseBD_est = zeros(K*Ns, Nr);
        for k = 1:K
            WphaseBD_est((k-1)*Ns + (1:Ns), :) = WBD_est((k-1)*Ns + (1:Ns), :) * (Wff_est(:, (k-1)*Lr + (1:Lr))');
        end

        TBD_est = TBD_est*sqrt(K*Ns)/norm(Trf_est*TBD_est, 'fro');
        TphaseBD_est = Trf_est*TBD_est;
        
        ratephaseBD_dft_err(isnr) = ratephaseBD_dft_err(isnr) + calRateMU(H, TphaseBD_est, WphaseBD_est, Rs, Ns*ones(K,1));    
        
        
    end
    fprintf('SNR_err = %d dB\n', SNR_err(isnr));
end

rateBD = rateBD/channNum;
ratephaseBD_dft = ratephaseBD_dft/channNum;
rateBD_err = rateBD_err/channNum;
ratephaseBD_dft_err = ratephaseBD_dft_err/channNum;

RAY_RATE_SET = [rateBD, rateBD_err, ratephaseBD_dft, ratephaseBD_dft_err];

%%%%% Plotting
figure
lw = 1.5;
ms = 5;
hold on
plot(SNR_err, abs(rateBD), 'k--', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR_err, abs(rateBD_err), 'k*', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR_err, abs(ratephaseBD_dft), 'r-', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR_err, abs(ratephaseBD_dft_err), 'ro', 'LineWidth', lw, 'MarkerSize', ms)
% plot(SNR, abs(ratephaseBD_svd_dft), 'ro', 'LineWidth', lw, 'MarkerSize', ms)
hold off

legend('Full-complexity BD', 'Full-complexity BD(err)', 'Hybrid BD', 'Hybrid BD(err)')
xlabel('Normalized channel estimation error power (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')
title(sprintf('Nt = %d, Nr = %d, K = %d, Ns = %d, Lr = %d, Lt = %d', ...
    Nt, Nr, K, Ns, Lr, Lt))
grid

filename = sprintf('BD-Nt%d-K%d-Nr%d-Ns%d-Lr%d-Lt%d', Nt, K, Nr, Ns, ...
     Lr,Lt);

save(filename, 'RAY_RATE_SET')

saveas(gcf, sprintf('BD-Nt%d-K%d-Nr%d-Ns%d-Lr%d-Lt%d', Nt, K, Nr, Ns, ...
    Lr,Lt)); 

toc
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

K = 8;      % User number
Ns = 2;     % #streams per user
Lr = 2;     % #chains per user
Lt = Lr*K;  % #chain at BS

% randn('state', 3);
channNum = 100;  % For each user
display('Channel Loading ...')
genH = zeros(Nr, Nt, channNum*K);
for ch = 1:channNum*K
    genH(:, :, ch) = (randn(Nr, Nt) + j*randn(Nr, Nt)) / sqrt(2);
end
display('Successful!');

SNR = -40 : 5 : 0;
nSNR = length(SNR);

j = sqrt(-1);

rateBD = zeros(nSNR, 1);        % full-complexity BD
ratephaseBD_dft = zeros(nSNR, 1);   % Phase + BB-BD


for isnr = 1 : nSNR
    P = 10^(SNR(isnr)/10);
    for ichannel = 1 : channNum   
        
        % ============= channel generation =====================
        H = zeros(K*Nr, Nt);
        Hcell = cell(K, 1);
        for k = 1:K
            Hcell{k} = genH(:, :, ichannel + channNum * (k-1));
            H((k-1)*Nr + (1:Nr), 1:Nt) = Hcell{k};
        end
       
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
        
    end
    fprintf('SNR = %d dB\n', SNR(isnr));
end

rateBD = rateBD/channNum;
ratephaseBD_dft = ratephaseBD_dft/channNum;

RAY_RATE_SET = [rateBD, ratephaseBD_dft];

%%%%% Plotting
figure
lw = 1.5;
ms = 5;
hold on
plot(SNR, abs(rateBD), 'k-*', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(ratephaseBD_dft), 'r-o', 'LineWidth', lw, 'MarkerSize', ms)
% plot(SNR, abs(ratephaseBD_svd_dft), 'ro', 'LineWidth', lw, 'MarkerSize', ms)
hold off

legend('Full-complexity BD', 'Hybrid BD')
xlabel('SNR (dB)')
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
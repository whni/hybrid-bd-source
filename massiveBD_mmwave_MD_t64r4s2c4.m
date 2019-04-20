% Study capacity performance of BD in massive MIMO
% Hybrid precoding VS full-complexity BD
% By Weiheng Ni, UVic, Oct. 27, 2013
% Modified July 16, 2014

tic; clear all; clc;

% =============================================
% ============ Parameter settings =============
% =============================================
Nt = [64, 1];
Nr = [4, 1];
type = 1;

if type == 1
    Nt = prod(Nt);
    Nr = prod(Nr);
end

K = 8;      % User number
Ns = 2;     % #streams per user
Lr = 4;     % #chains per user
Lt = Lr*K;  % #chain at BS

% randn('state', 3);
Ncls = 8; 
Nray = 10;
channNum = 2;  % For each user
display('Channel Loading ...')
genH = zeros(prod(Nr), prod(Nt), channNum*K);
genAlpha = zeros(Ncls*Nray, channNum*K);
genAt = zeros(prod(Nt), Ncls*Nray, channNum*K);
genAr = zeros(prod(Nr), Ncls*Nray, channNum*K);
for k = 1:K
    if type == 1
        DirTX = [-60 60]*pi/180 + 2*pi*rand();
        DirRX = [-180 180]*pi/180;
    else
        DirTX = [-60 60; -60 60]*pi/180 + 2*pi*rand();
        DirRX = [-180 180; -180 180]*pi/180;
    end
    [genH(:, :, (1:channNum)+(k-1)*channNum), genAlpha(:, (1:channNum)+(k-1)*channNum), ...
        genAt(:, :, (1:channNum)+(k-1)*channNum), genAr(:,:,(1:channNum)+(k-1)*channNum)]...
        = channelSet(Nt, Nr, Ncls, Nray, channNum, DirTX, DirRX);
end
display('Successful!');

Nt = prod(Nt);
Nr = prod(Nr);

SNR = -40 : 5 : 0;
nSNR = length(SNR);

j = sqrt(-1);

rateBD = zeros(nSNR, 1);        % full-complexity BD
rateSpa = zeros(nSNR, 1);       % Hybrid BD based on sparse precoding
rateRxMmse = zeros(nSNR, 1);    % MMSE at RX
ratephaseBD_ar = zeros(nSNR, 1);   % Phase + BB-BD
ratephaseBD_dft = zeros(nSNR, 1);   % Phase + BB-BD
ratephaseBD_svd_dft = zeros(nSNR, 1);
ratephaseBD_svd_ar = zeros(nSNR, 1);
rateBD_MD = zeros(nSNR, 1);

for isnr = 1 : nSNR
    P = 10^(SNR(isnr)/10);
    for ichannel = 1 : channNum   
        
        % ============= channel generation =====================
        H = zeros(K*Nr, Nt);
        Hcell = cell(K, 1);
        At = cell(K, 1);
        Ar = cell(K, 1);
        Alpha = cell(K, 1);
        for k = 1:K
            Hcell{k} = genH(:, :, ichannel + channNum * (k-1));
            H((k-1)*Nr + (1:Nr), 1:Nt) = Hcell{k};
            At{k} = genAt(:, :, ichannel + channNum * (k-1));
            Ar{k} = genAr(:, :, ichannel + channNum * (k-1));
            Alpha{k} = genAlpha(:, ichannel + channNum * (k-1));
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

        
        % ================= BD + MD precoding =====================
        [TBD, WBD] = CalPrecoderBD(H, K, Ns);% BD precoding & combining
        [Trf, Tbb] = general_decomp(TBD, K*Ns, Lt);
        HH = H*TBD;
        Wff = cell(K, 1);
        Wbb = cell(K, 1);
        %WW = zeros(K*Ns, K*Nr);
        for k = 1:K
            WW((k-1)*Ns + (1:Ns), (k-1)*Nr + (1:Nr)) = WBD((k-1)*Ns + (1:Ns), 1:Nr);
            [Wff{k}, Wbb{k}] = general_decomp(WBD((k-1)*Ns + (1:Ns), 1:Nr)', Ns, Lr);
        end
        
        vec_power_gain = diag(abs(WW*HH)).^2; 
        vec_base_level = ones(K*Ns, 1) ./ vec_power_gain;
        vec_power_alloc = water_filling(P, vec_base_level);
        Fwf = diag(sqrt(K*Ns*vec_power_alloc/P));
        Tbb = Tbb*sqrt(K*Ns)/norm(Trf*Tbb*Fwf, 'fro');
        TphaseBD = Trf*Tbb*Fwf;
        
        WphaseBD = zeros(K*Ns, Nr);
        for k = 1:K
            WphaseBD((k-1)*Ns + (1:Ns), :) = (Wff{k} * Wbb{k})';
        end

        rateBD_MD(isnr) = rateBD_MD(isnr) + calRateMU(H, TphaseBD, WphaseBD, Rs, Ns*ones(K,1));
        
        
        
        % ================= Sparse Hybrid BD ================
        Tf = cell(K, 1);
        Tb = cell(K, 1);
        Tspa = zeros(Nt, K*Ns);
        Wf = cell(K, 1);
        Wb = cell(K, 1);
        Wspa = zeros(K*Ns, Nr);
        for k = 1:K
            [Tf{k}, Tb{k}] = CalSparsePrecoder(At{k}, TBD(:, (k-1)*Ns + (1:Ns)), Lt/K); 
            Tspa(:, (k-1)*Ns + (1:Ns)) = Tf{k}*Tb{k};
            [Wf{k}, Wb{k}] = CalSparsePrecoder(Ar{k}, WBD((k-1)*Ns + (1:Ns), :)', Lr);
            Wspa((k-1)*Ns + (1:Ns), :) = (Wf{k}*Wb{k})';
        end
        
        rateSpa(isnr) = rateSpa(isnr) + calRateMU(H, Tspa, Wspa, Rs, Ns*ones(K,1));

        
        % ================ BD at TX, MMSE at RX ===============
        Tmmse = Tspa;
        Wmmse = cell(K, 1);
        Eyy = cell(K, 1);
        WRxMmse = zeros(K*Ns, Nr);
        for k = 1:K
            zero_eye = zeros(Ns, K*Ns);
            zero_eye(:, (k-1)*Ns + (1:Ns)) = eye(Ns);
            Wmmse{k} = 1/sqrt(P) * zero_eye * inv(Tmmse'*Hcell{k}'*Hcell{k}*Tmmse + K*Ns/P*eye(K*Ns)) *Tmmse' *Hcell{k}';
            Wmmse{k} = Wmmse{k}';
            Eyy{k} = (P/(K*Ns)) * Hcell{k} * Tmmse*Tmmse' * Hcell{k}' + eye(Nr);
            [Wf{k}, Wb{k}] = calMmseCombiner(Ar{k}, Wmmse{k}, Eyy{k}, Lr);
            WRxMmse((k-1)*Ns + (1:Ns), :) = (Wf{k}*Wb{k})';
        end

        rateRxMmse(isnr) = rateRxMmse(isnr) + calRateMU(H, Tmmse, WRxMmse, Rs, Ns*ones(K,1));      
        
        
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
rateSpa = rateSpa/channNum;
rateRxMmse = rateRxMmse/channNum;
ratephaseBD_ar = ratephaseBD_ar/channNum;
ratephaseBD_dft = ratephaseBD_dft/channNum;
ratephaseBD_svd_dft = ratephaseBD_svd_dft/channNum;
ratephaseBD_svd_ar = ratephaseBD_svd_ar/channNum;
rateBD_MD = rateBD_MD/channNum;
%%%%% Plotting

rateMat = [rateBD, rateRxMmse, ratephaseBD_dft, rateBD_MD];

figure
lw = 1.5;
ms = 5;
hold on
plot(SNR, abs(rateBD), 'k--', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(ratephaseBD_dft), 'r-', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rateRxMmse), 'b-.', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rateBD_MD), 'g*-', 'LineWidth', lw, 'MarkerSize', ms)

hold off
legend('Full-Complexity BD', 'Hybrid BD', 'Spatially Sparse Coding')
xlabel('SNR (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')
title(sprintf('Nt = %d, Nr = %d, K = %d, Ns = %d, Lr = %d, Lt = %d, Ncls = %d, Nray = %d', ...
    Nt, Nr, K, Ns, Lr, Lt, Ncls, Nray))
grid

filename = sprintf('mmwave-MD-Nt%d-K%d-Nr%d-Ns%d-Lr%d-Lt%d', Nt, K, Nr, Ns, Lr, Lt);

save(filename, 'rateMat')

% saveas(gcf, sprintf('BD-Nt%d-K%d-Nr%d-Ns%d-Lr%d-Lt%d-Ncls%d-Nray%d', Nt, K, Nr, Ns, ...
%     Lr,Lt,Ncls,Nray)); 

toc
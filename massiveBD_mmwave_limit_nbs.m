% Study capacity performance of BD in massive MIMO
% Hybrid precoding VS full-complexity BD
% By Weiheng Ni, UVic, Oct. 27, 2013
% Modified July 16, 2014

tic; clear all; clc;

% =============================================
% ============ Parameter settings =============
% =============================================

NtSet = 32:64:512;
nNt = length(NtSet);

rateBD = zeros(nNt, 1);        % full-complexity BD
rateSpa = zeros(nNt, 1);       % Hybrid BD based on sparse precoding
rateRxMmse = zeros(nNt, 1);    % MMSE at RX
ratephaseBD_ar = zeros(nNt, 1);   % Phase + BB-BD
ratephaseBD_dft = zeros(nNt, 1);   % Phase + BB-BD
ratephaseBD_svd_dft = zeros(nNt, 1);
ratephaseBD_svd_ar = zeros(nNt, 1);
rateLF = zeros(nNt, 1);
rateAL = zeros(nNt, 1);

Nr = 64;
type = 1;

K = 2;      % User number
Ns = 1;     % #streams per user
Lr = 1;     % #chains per user
Lt = Lr*K;  % #chain at BS

% randn('state', 3);
SNR = 0; 
P = 10^(SNR/10);

Ncls = 1; 
Nray = 1;
channNum = 200;  % For each user
large_fa_loss = 0.2; %  large scale fading lower bound

for iNt = 1:1:nNt
    Nt = NtSet(iNt);
    Nr = Nt / (2*K);
    display('Channel Loading ...')
    genH = zeros(Nr, Nt, channNum*K);
    genAlpha = zeros(Ncls*Nray, channNum*K);
    genAt = zeros(Nt, Ncls*Nray, channNum*K);
    genAr = zeros(Nr, Ncls*Nray, channNum*K);
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
            = channelSet(Nt, Nr, Ncls, Nray, channNum, DirTX, DirRX, large_fa_loss);
    end
    display('Successful!');
    
    j = sqrt(-1);

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
        rateBD(iNt) = rateBD(iNt) + calRateMU(H, TBDwf, WBD, Rs, Ns*ones(K,1));
        
        
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
        
        rateSpa(iNt) = rateSpa(iNt) + calRateMU(H, Tspa, Wspa, Rs, Ns*ones(K,1));
        
        % ================= Limited Feedback ================
        Wff = zeros(Nr, K*Lr);
        Trf = zeros(Nt, K*Lr);
        Heq = zeros(K*Lr, Nt);
        Heq2 = zeros(K*Lr, Nt);
        
        for k = 1:K
            [Trf(:, (k-1)*Lr + (1:Lr)), Wff(:, (k-1)*Lr + (1:Lr))] = limited_feedback_bases(Hcell{k}, 1, At{k}, Ar{k});
            Heq((k-1)*Lr + (1:Lr), :) = Wff(:, (k-1)*Lr + (1:Lr))' * Hcell{k} ;
        end
        
        Hbb = Heq*Trf;
        TBB = Hbb' / (Hbb * Hbb');
        
        Wlf = zeros(K*Ns, Nr);
        for k = 1:K
            Wlf((k-1)*Ns + (1:Ns), :) = (Wff(:, (k-1)*Lr + (1:Lr))');
        end
        
        vec_power_gain = diag(abs(Hbb*TBB)).^2;
        vec_base_level = ones(K*Ns, 1) ./ vec_power_gain;
        vec_power_alloc = water_filling(P, vec_base_level);
        Fwf = diag(sqrt(K*Ns*vec_power_alloc/P));
        
        Tlf = Trf*TBB*Fwf*sqrt(K*Ns)/norm(Trf*TBB*Fwf, 'fro');
        
        rateLF(iNt) = rateLF(iNt) + calRateMU(H, Tlf, Wlf, Rs, Ns*ones(K,1));
        
        
        % ==================phase + BD (DFT)====================
        Wff = zeros(Nr, K*Lr);
        Wff2 = zeros(Nr, K*Lr);
        
        Heq = zeros(K*Lr, Nt);
        for k = 1:K
            Wff(:, (k-1)*Lr + (1:Lr)) = select_wf_bases(Hcell{k}, 1, 1/sqrt(Nr) * dftmtx(Nr), Lr);
            Heq((k-1)*Lr + (1:Lr), :) = Wff(:, (k-1)*Lr + (1:Lr))' * Hcell{k};
        end
%         Wff = [Ar{:}];
%         for k = 1:K
%             Heq((k-1)*Lr + (1:Lr), :) = Wff(:, (k-1)*Lr + (1:Lr))' * Hcell{k};
%         end
        
        Wff2 = [Ar{:}];
        for k = 1:K
            Heq2((k-1)*Lr + (1:Lr), :) = Wff2(:, (k-1)*Lr + (1:Lr))' * Hcell{k};
        end
        Trf = exp(-j*angle(Heq2)') / sqrt(Nt);
        Hbb2 = Heq2 * Trf;
        
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
        
        ratephaseBD_dft(iNt) = ratephaseBD_dft(iNt) + calRateMU(H, TphaseBD, WphaseBD, Rs, Ns*ones(K,1));
        
        
        %========================alpha=====================
        Wal = ones(K*Lr, 1);
        %Aal = diag(diag(Hbb2));
        Aal = diag([Alpha{:}]);
        
        vec_power_gain = diag(abs(Aal)).^2;
        vec_base_level = ones(K*Ns, 1) ./ vec_power_gain;
        vec_power_alloc = water_filling(P, vec_base_level);
        Fwf = diag(sqrt(K*Ns*vec_power_alloc/P));
        Tal = Fwf;
        
        rateAL(iNt) = rateAL(iNt) + calRateMU(Aal, Tal, Wal, Rs, Ns*ones(K, 1));
        
        
    end
    fprintf('Nt = %d\n', NtSet(iNt));
end

rateBD = rateBD/channNum;
rateSpa = rateSpa/channNum;
rateRxMmse = rateRxMmse/channNum;
ratephaseBD_ar = ratephaseBD_ar/channNum;
ratephaseBD_dft = ratephaseBD_dft/channNum;
ratephaseBD_svd_dft = ratephaseBD_svd_dft/channNum;
ratephaseBD_svd_ar = ratephaseBD_svd_ar/channNum;
rateLF = rateLF/channNum;
 rateAL =  rateAL/channNum;
%%%%% Plotting

rateMat = [rateBD, rateSpa, ratephaseBD_dft, rateLF];

figure
lw = 1.5;
ms = 5;
hold on
plot(NtSet, abs(rateBD), 'k-*', 'LineWidth', lw, 'MarkerSize', ms)
plot(NtSet, abs(ratephaseBD_dft), 'b--', 'LineWidth', lw, 'MarkerSize', ms)
plot(NtSet, abs(rateAL), 'bv', 'LineWidth', lw, 'MarkerSize', ms)
plot(NtSet, abs(rateLF), 'g-o', 'LineWidth', lw, 'MarkerSize', ms)

hold off
legend('Full-Complexity BD', 'Hybrid BD', 'Upper Bound', 'Limited Feedback')
xlabel('N_{BS}')
ylabel('Sum spectral efficiency (bps/Hz)')
title(sprintf('Nt = %d, Nr = %d, K = %d, Ns = %d, Lr = %d, Lt = %d, Ncls = %d, Nray = %d', ...
    Nt, Nr, K, Ns, Lr, Lt, Ncls, Nray))
grid

filename = sprintf('ALpha-Nt%d-K%d-Nr%d-Ns%d-Lr%d-Lt%d-Ncls%d-Nray%d', Nt, K, Nr, Ns, ...
     Lr,Lt,Ncls,Nray);

save(filename, 'rateMat')

% saveas(gcf, sprintf('BD-Nt%d-K%d-Nr%d-Ns%d-Lr%d-Lt%d-Ncls%d-Nray%d', Nt, K, Nr, Ns, ...
%     Lr,Lt,Ncls,Nray)); 

toc
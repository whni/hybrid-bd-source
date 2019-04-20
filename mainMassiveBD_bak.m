% Study capacity performance of BD in massive MIMO
% Hybrid precoding VS full-complexity BD
% By Weiheng Ni, UVic, Oct. 27, 2013
% Modified July 16, 2014

tic; clear all; clc;

% ============ Parameter settings =============
% =============================================
Nt = 128;
Nr = 32;
K = 3;      % User number
Ns = 2;     % #streams per user
Lr = 4;     % #chains per user
Lt = Lr*K;  % #chain at BS

% randn('state', 3);
Ncls = 8; 
Nray = 10;
channNum = 10;  % For each user
display('Channel Loading ...')
% load E:\MatlabData\DIRchannel-Nt64-Nr16-Ncls8-Nray10-channNum1000.mat; Ncls = 8; Nray = 10;
% load E:\MatlabData\OMNIchannel-Nt64-Nr16-Ncls8-Nray10-channNum1000.mat; Ncls = 8; Nray = 10;
% load ChannelData\OMNIchannel-Nt128-Nr16-Ncls10-Nray8-channNum1000.mat; 
[genH, genAlpha, genAt, genAr] = channelSet(Nt, Nr, Ncls, Nray, channNum*K);
display('Successful!')


% Nb = Ncls*Nray; % #bases for each user
% at_bases = zeros(Nt, K*Nb); % bases for all users

SNR = -40 : 5 : 0;
nSNR = length(SNR);

j = sqrt(-1);

% ============== channelSet to get a bunch of channels =================
% [genH, genAlpha, genAt, genAr] = channelSet(sqrt(Nt)*ones(2,1), ...
%     sqrt(Nr)*ones(2,1), Ncls, Nray, 2*channNum);
% ======================================================================

rateBD = zeros(nSNR, 1);        % full-complexity BD
rateSpa = zeros(nSNR, 1);       % Hybrid BD based on sparse precoding
rateRxMmse = zeros(nSNR, 1);    % MMSE at RX
ratephaseBD = zeros(nSNR, 1);   % Phase + BB-BD

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
        
%         H1 = genH(:, :, ichannel);% channel for UE1
%         At1 = genAt(:, :, ichannel);
%         Ar1 = genAr(:, :, ichannel);
%         Alpha1 = genAlpha(:, ichannel);
%         
%         H2 = genH(:, :, ichannel+channNum);% channel for UE2
%         At2 = genAt(:, :, ichannel+channNum);
%         Ar2 = genAr(:, :, ichannel+channNum);    
%         Alpha2 = genAlpha(:, ichannel+channNum);
%         
%         H3 = genH(:, :, ichannel+2*channNum);% channel for UE3
%         At3 = genAt(:, :, ichannel+2*channNum);
%         Ar3 = genAr(:, :, ichannel+2*channNum);
%         Alpha3 = genAlpha(:, ichannel+2*channNum);
%         
%         H4 = genH(:, :, ichannel+3*channNum);% channel for UE4
%         At4 = genAt(:, :, ichannel+3*channNum);
%         Ar4 = genAr(:, :, ichannel+3*channNum);    
%         Alpha4 = genAlpha(:, ichannel+3*channNum);
%         
%         H = [H1; H2; H3; H4];
        
%         at_bases = [At1, At2, At3, At4];
%         at_angle_set = gen_bases(H, K, Nb);
%         for kb = 1:(K*Nb)
%             at_bases(:, kb) = exp(j*2*pi*0.5*sin(at_angle_set(kb)) *([0:(Nt-1)]') )/sqrt(Nt);
%         end
%         at_bases = select_at_bases(H, K, [At1, At2]);


        Rs = P/(K*Ns)*eye(K*Ns);
      
        % ================= BD precoding =====================
        [TBD, WBD] = CalPrecoderBD(H, K, Ns);% BD precoding & combining
        HH = H*TBD;
        WW = zeros(K*Ns, K*Nr);
        for k = 1:K
            WW((k-1)*Ns + (1:Ns), (k-1)*Nr + (1:Nr)) = WBD((k-1)*Ns + (1:Ns), 1:Nr);
        end
%         WW = [WBD(1:Ns, :), zeros(2, 48); zeros(2, 16),WBD(Ns+1:2*Ns, :), zeros(2, 32); ...
%             zeros(2, 32), WBD(1+2*Ns:3*Ns, :), zeros(2, 16); zeros(2, 48), WBD(3*Ns+1:4*Ns, :) ];
        vec_power_gain = diag(abs(WW*HH)).^2; 
        vec_base_level = ones(K*Ns, 1) ./ vec_power_gain;
        vec_power_alloc = water_filling(P, vec_base_level);
        Fwf = diag(sqrt(K*Ns*vec_power_alloc/P));
        TBDwf = TBD*Fwf;
        rateBD(isnr) = rateBD(isnr) + calRateMU(H, TBDwf, WBD, Rs, Ns*ones(K,1));


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
        
%         [Tf1, Tb1] = CalSparsePrecoder(At1, TBD(:, 1:Ns), Lt/K); 
%         [Tf2, Tb2] = CalSparsePrecoder(At2, TBD(:, (Ns+1):(2*Ns)), Lt/K);
%         [Tf3, Tb3] = CalSparsePrecoder(At3, TBD(:, (2*Ns+1):(3*Ns)), Lt/K); 
%         [Tf4, Tb4] = CalSparsePrecoder(At4, TBD(:, (3*Ns+1):(4*Ns)), Lt/K); 
% %         Tf = [Tf1, Tf2, Tf3, Tf4];
% %         Tb = [Tb1 zeros(Lt/2, Ns); zeros(Lt/2, Ns), Tb2];
%         Tspa = [Tf1*Tb1, Tf2*Tb2, Tf3*Tb3, Tf4*Tb4];% TX precoder
%         [Wf1, Wb1] = CalSparsePrecoder(Ar1, WBD(1:Ns,:)', Lr);
%         [Wf2, Wb2] = CalSparsePrecoder(Ar2, WBD((1+Ns):2*Ns,:)', Lr);
%         [Wf3, Wb3] = CalSparsePrecoder(Ar3, WBD((1+2*Ns):3*Ns,:)', Lr);
%         [Wf4, Wb4] = CalSparsePrecoder(Ar4, WBD((1+3*Ns):4*Ns,:)', Lr);
%         W1 = (Wf1*Wb1)';
%         W2 = (Wf2*Wb2)';
%         W3 = (Wf3*Wb3)';
%         W4 = (Wf4*Wb4)';
%         Wspa = [W1; W2; W3; W4];% RX combiner
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
%         Wmmse1 = 1/sqrt(P)*[eye(Ns) zeros(Ns) zeros(Ns) zeros(Ns)]*inv(Tmmse'*H1'*H1*Tmmse + K*Ns/P*eye(K*Ns))*Tmmse'*H1';
%         Wmmse2 = 1/sqrt(P)*[zeros(Ns) eye(Ns) zeros(Ns) zeros(Ns)]*inv(Tmmse'*H2'*H2*Tmmse + K*Ns/P*eye(K*Ns))*Tmmse'*H2';
%         Wmmse3 = 1/sqrt(P)*[zeros(Ns) zeros(Ns) eye(Ns) zeros(Ns)]*inv(Tmmse'*H3'*H3*Tmmse + K*Ns/P*eye(K*Ns))*Tmmse'*H3';
%         Wmmse4 = 1/sqrt(P)*[zeros(Ns) zeros(Ns) zeros(Ns) eye(Ns)]*inv(Tmmse'*H4'*H4*Tmmse + K*Ns/P*eye(K*Ns))*Tmmse'*H4';
%         Wmmse1 = Wmmse1'; Wmmse2 = Wmmse2';Wmmse3 = Wmmse3'; Wmmse4 = Wmmse4';
% 
%         Eyy1 = (P/(K*Ns))*H1*Tmmse*Tmmse'*H1' + eye(Nr);
%         Eyy2 = (P/(K*Ns))*H2*Tmmse*Tmmse'*H2' + eye(Nr);
%         Eyy3 = (P/(K*Ns))*H3*Tmmse*Tmmse'*H3' + eye(Nr);
%         Eyy4 = (P/(K*Ns))*H4*Tmmse*Tmmse'*H4' + eye(Nr);
%         [Wf1, Wb1] = calMmseCombiner(Ar1, Wmmse1, Eyy1, Lr);
%         [Wf2, Wb2] = calMmseCombiner(Ar2, Wmmse2, Eyy2, Lr);
%         [Wf3, Wb3] = calMmseCombiner(Ar3, Wmmse3, Eyy3, Lr);
%         [Wf4, Wb4] = calMmseCombiner(Ar4, Wmmse4, Eyy4, Lr);
%         WRxMmse = [(Wf1*Wb1)'; (Wf2*Wb2)'; (Wf3*Wb3)'; (Wf4*Wb4)'];
        rateRxMmse(isnr) = rateRxMmse(isnr) + calRateMU(H, Tmmse, WRxMmse, Rs, Ns*ones(K,1));      
        
        
        % ==================phase + BD ====================
%         [~, max_alpha_pos1] = sort(Alpha1.^2, 'descend');
%         Wf1 = Ar1(:, max_alpha_pos1(1:Lr));
%         [~, max_alpha_pos2] = sort(Alpha2.^2, 'descend');
%         Wf2 = Ar2(:, max_alpha_pos2(1:Lr));
%         Wf1 = Wff1;
%         Wf2 = Wff2;
%         Wf3 = Wff3;
%         Wf4 = Wff4;
        Ar_set = zeros(Nr, K*Ncls*Nray);
        for k = 1:K;
            Ar_set(:, (k-1)*Ncls*Nray + (1:Ncls*Nray)) = Ar{k};
        end
        Wff = select_wf_bases(H, K, Ar_set, Lr);
%         Wf = Wff;
%         Wf1 = Wff(:, 1:Lr);
%         Wf2 = Wff(:, Lr+1:2*Lr);
%         Wf3 = Wff(:, 2*Lr+1:3*Lr);
%         Wf4 = Wff(:, 3*Lr+1:4*Lr);

        Heq = zeros(K*Lr, Nt);
        for k = 1:K
            Heq((k-1)*Lr + (1:Lr), :) = Wff(:, (k-1)*Lr + (1:Lr))' * Hcell{k};
        end
%         Heq1 = Wf1' * H1;
%         Heq2 = Wf2' * H2;
%         Heq3 = Wf3' * H3;
%         Heq4 = Wf4' * H4;
%         Heq = [Heq1; Heq2; Heq3; Heq4];
        Trf = exp(-j*angle(Heq)') / sqrt(Nt);
        Hbb = Heq * Trf;
        [TBD, WBD] = CalPrecoderBD(Hbb, K, Ns);
        
        HH = Hbb*TBD;
        WW = zeros(K*Ns, K*Lr);
        for k = 1:K
            WW((k-1)*Ns + (1:Ns), (k-1)*Lr + (1:Lr)) = WBD((k-1)*Ns + (1:Ns), :);
        end
%         WW = [WBD(1:Ns, :), zeros(2, 12); zeros(2, 4),WBD(Ns+1:2*Ns, :), zeros(2, 8); ...
%             zeros(2, 8), WBD(1+2*Ns:3*Ns, :), zeros(2, 4); zeros(2, 12), WBD(3*Ns+1:4*Ns, :) ];
        vec_power_gain = diag(abs(WW*HH)).^2; 
        vec_base_level = ones(K*Ns, 1) ./ vec_power_gain;
        vec_power_alloc = water_filling(P, vec_base_level);
        Fwf = diag(sqrt(K*Ns*vec_power_alloc/P));
        
        WphaseBD = zeros(K*Ns, Nr);
        for k = 1:K
            WphaseBD((k-1)*Ns + (1:Ns), :) = WBD((k-1)*Ns + (1:Ns), :) * (Wff(:, (k-1)*Lr + (1:Lr))');
        end
%         WphaseBD = [WBD(1:Ns, :)*(Wf1'); WBD(Ns+1:2*Ns, :)*(Wf2'); WBD(1+2*Ns:3*Ns, :)*(Wf3'); WBD(3*Ns+1:4*Ns, :)*(Wf4')];
        
        TBD = TBD*sqrt(K*Ns)/norm(Trf*TBD, 'fro');
        TphaseBD = Trf*TBD*Fwf;
        
        ratephaseBD(isnr) = ratephaseBD(isnr) + calRateMU(H, TphaseBD, WphaseBD, Rs, Ns*ones(K,1));    
        
    end
    fprintf('SNR = %d dB\n', SNR(isnr));
end

rateBD = rateBD/channNum;
rateSpa = rateSpa/channNum;
rateRxMmse = rateRxMmse/channNum;
ratephaseBD = ratephaseBD/channNum;

%%%%% Plotting
figure
lw = 1.5;
ms = 5;
hold on
plot(SNR, abs(rateBD), 'k-*', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rateSpa), 'b-o', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(rateRxMmse), 'r-^', 'LineWidth', lw, 'MarkerSize', ms)
plot(SNR, abs(ratephaseBD), 'g--', 'LineWidth', lw, 'MarkerSize', ms)
hold off
legend('BD', 'Spatially Sparse Hybrid', 'MMSE at RX', 'phase BD')
xlabel('SNR (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')
title(sprintf('Nt = %d, Nr = %d, K = %d, Ns = %d, Lr = %d, Lt = %d, Ncls = %d, Nray = %d', ...
    Nt, Nr, K, Ns, Lr, Lt, Ncls, Nray))
grid
% saveas(gcf, sprintf('BD-Nt%d-K%d-Nr%d-Ns%d-Lr%d-Lt%d-Ncls%d-Nray%d', Nt, K, Nr, Ns, ...
%     Lr,Lt,Ncls,Nray)); 

toc
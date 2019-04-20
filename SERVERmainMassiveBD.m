% Study capacity performance of BD in massive MIMO
% Hybrid precoding VS full-complexity BD
% By Le Nsiang, UVic, Oct. 27, 2013
% Modified July 16, 2014

tic; clear all; clc;
% randn('state', 3);

display('Channel Loading ...')
% load OMNIchannel-Nt64-Nr16-Ncls8-Nray10-channNum1000.mat; Ncls = 8; Nray = 10;
% load OMNIchannel-Nt64-Nr16-Ncls10-Nray1-channNum1000.mat; Ncls = 10; Nray = 1;
load OMNIchannel-Nt256-Nr64-Ncls10-Nray1-channNum1000.mat; Ncls = 10; Nray = 1;
display('Successful!')

% ============ Parameter settings =============
% =============================================
Nt = 256; % TX antenna number, MU-MIMO
K = 2; % User number
Nr = 64; % RX antenna number
Ns = 2;% #streams per user
Lr = 2;% #chains per user
Lt = 6;% #chain at BS


SNR = -30 : 5 : 0;
nSNR = length(SNR);
channNum = 5e2; % channel realization times

% ============== channelSet to get a bunch of channels =================
% [genH, genAlpha, genAt, genAr] = channelSet(sqrt(Nt)*ones(2,1), ...
%     sqrt(Nr)*ones(2,1), Ncls, Nray, 2*channNum);
% ======================================================================

rateBD = zeros(nSNR, 1);
rateSpa = zeros(nSNR, 1);% Hybrid BD based on sparse precoding
rateRxMmse = zeros(nSNR, 1);% MMSE at RX

for isnr = 1 : nSNR
    P = 10^(SNR(isnr)/10);
    for ichannel = 1 : channNum   
        
        % ============= channel generation =====================    
        H1 = genH(:, :, ichannel);% channel for UE1
        At1 = genAt(:, :, ichannel);
        Ar1 = genAr(:, :, ichannel);
        
        H2 = genH(:, :, ichannel+channNum);% channel for UE2
        At2 = genAt(:, :, ichannel+channNum);
        Ar2 = genAr(:, :, ichannel+channNum);    
        H = [H1; H2];
    
        Rs = P/(K*Ns)*eye(K*Ns);
      
        % ================= BD precoding =====================
        [TBD, WBD] = CalPrecoderBD(H, K, Ns);% BD precoding & combining
        rateBD(isnr) = rateBD(isnr) + calRateMU(H, TBD, WBD, Rs, Ns*ones(K,1));

        % ================= Sparse Hybrid BD ================
        [Tf Tb] = CalSparsePrecoder([At1 At2], TBD, Lt);
        Tspa = Tf*Tb;% TX precoder
        [Wf1, Wb1] = CalSparsePrecoder(Ar1, WBD(1:Ns,:)', Lr);
        [Wf2, Wb2] = CalSparsePrecoder(Ar2, WBD((1+Ns):2*Ns,:)', Lr);
        W1 = (Wf1*Wb1)';
        W2 = (Wf2*Wb2)';
        Wspa = [W1; W2];% RX combiner
        rateSpa(isnr) = rateSpa(isnr) + calRateMU(H, Tspa, Wspa, Rs, Ns*ones(K,1));
        
        % ================ BD at TX, MMSE at RX ===============
        Tmmse = Tspa;
        Wmmse1 = 1/sqrt(P)*[eye(Ns) zeros(Ns)]*inv(Tmmse'*H1'*H1*Tmmse + 2*Ns/P*eye(2*Ns))*Tmmse'*H2';
        Wmmse2 = 1/sqrt(P)*[zeros(Ns) eye(Ns)]*inv(Tmmse'*H2'*H2*Tmmse + 2*Ns/P*eye(2*Ns))*Tmmse'*H2';
        Wmmse1 = Wmmse1'; Wmmse2 = Wmmse2';
        Eyy1 = 1/(2*Ns)*H1*Tmmse*Tmmse'*H1';
        Eyy2 = 1/(2*Ns)*H2*Tmmse*Tmmse'*H2';
        [Wf1, Wb1] = calMmseCombiner(Ar1, Wmmse1, Eyy1, Lr);
        [Wf2, Wb2] = calMmseCombiner(Ar2, Wmmse2, Eyy2, Lr);
        WRxMmse = [(Wf1*Wb1)'; (Wf2*Wb2)'];
        rateRxMmse(isnr) = rateRxMmse(isnr) + calRateMU(H, Tmmse, WRxMmse, Rs, Ns*ones(K,1));      
    end
    isnr
end

rateBD = rateBD/channNum;
rateSpa = rateSpa/channNum;
rateRxMmse = rateRxMmse/channNum;

%%%%% Plotting
figure
lw = 1.5;
ms = 6;
plot(SNR, abs(rateBD), 'k-*', 'LineWidth', lw, 'MarkerSize', ms)
hold on
plot(SNR, abs(rateSpa), 'b-o', 'LineWidth', lw, 'MarkerSize', ms)
hold on
plot(SNR, abs(rateRxMmse), 'r-^', 'LineWidth', lw, 'MarkerSize', ms)
hold off
legend('BD', 'Spatially Sparse Hybrid', 'MMSE at RX')
xlabel('SNR (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')
title(sprintf('Nt = %d, Nr = %d, K = %d, Ns = %d, Lr = %d, Lt = %d, Ncls = %d, Nray = %d', ...
    Nt, Nr, K, Ns, Lr, Lt, Ncls, Nray))
grid
saveas(gcf, sprintf('BD-Nt%d-K%d-Nr%d-Ns%d-Lr%d-Lt%d-Ncls%d-Nray%d', Nt, K, Nr, Ns, ...
    Lr,Lt,Ncls,Nray)); 

toc
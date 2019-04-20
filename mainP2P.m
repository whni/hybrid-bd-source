% Massive P2P MIMO, hybrid RF-baseband precoding
% By Le Liang, UVic, July 2, 2014

clear all; clc; tic;
% randn('state', 5);

display('Loading channel ...')
% load E:\MatlabData\mmWaveChanH.mat % Leyuan's channel data
load E:\MatlabData\DIRchannel-Nt64-Nr16-Ncls8-Nray10-channNum1000.mat
% load E:\MatlabData\OMNIchannel-Nt64-Nr16-Ncls8-Nray10-channNum1000.mat
display('Successful!')

Nt = 64;
Nr = 16;
Ns = 2;
L = 2; % #chains at TX/RX

Ncls = 8;% #clusters  
Nray = 10;% #rays per cluster

SNR = -30 : 5 : 0;
nSNR = length(SNR);
channNum = 1e3;

% ============== channelSet to get a bunch of channels =================
% [genH, genAlpha, genAt, genAr] = channelSet(sqrt(Nt)*ones(2,1), ...
%     sqrt(Nr)*ones(2,1), Ncls, Nray, channNum);
% ======================================================================

rateSvdEq = zeros(nSNR, 1);
rateSpaEq = zeros(nSNR, 1);% spatially sparse precoding of Heath's paper
rateHybEq = zeros(nSNR, 1);
rateSvdWf = zeros(nSNR, 1);% water-filling capacity
rateSpaWf = zeros(nSNR, 1);
rateHybWf = zeros(nSNR, 1);

for isnr = 1 : nSNR
    P = 10^(SNR(isnr)/10);
    for ichannel = 1 : channNum
        
        H = genH(:, :, ichannel);% channel instantiation from channel set
        At = genAt(:, :, ichannel);
        Ar = genAr(:, :, ichannel);   
        
%         H = chanH(:, :, ichannel); % Leyuan's channel data
%         At = reshape(aAntArrayTX(:, :, :, ichannel), [Nt, Ncls*Nray]);
%         Ar = reshape(aAntArrayRX(:, :, :, ichannel), [Nr, Ncls*Nray]);
%         
        [U, S, V] = svd(H);
        
        %==================================================================
        % =========== optimum uncontrained precoding based on SVD =========
        %==================================================================
        Topt = V(:, 1:Ns);% TX precoding matrix
        Wopt = U(:, 1:Ns)';% RX combining matrix
        
        % ====== water-filling power allocation
%         gain = diag(S);
%         RsWf = waterfill(P, gain(1:Ns));% power allocation
%         RsWf = diag(RsWf);
        
        % ====== equal power allocation 
        RsEq = P/Ns*eye(Ns);
        
        Rn = Wopt*Wopt';
        WHT = Wopt*H*Topt;
        rateTmpEq = log2(det(eye(Ns) + inv(Rn)*WHT*RsEq*WHT'));
        rateSvdEq(isnr) = rateSvdEq(isnr) + rateTmpEq;
%         rateTmpWf = log2(det(eye(Ns) + inv(Rn)*WHT*RsWf*WHT'));
%         rateSvdWf(isnr) = rateSvdWf(isnr) + rateTmpWf;
        
        %==================================================================
        % ============ Spatially sparse precoding =========================
        %==================================================================
        [Tf, Tb] = CalSparsePrecoder(At, V(:, 1:Ns), L);% equal power
%         [TfWf, TbWf] = CalSparsePrecoder(At, V(:, 1:Ns)*sqrt(RsWf), L);% waterfillling
        [Wf, Wb] = CalSparsePrecoder(Ar, U(:, 1:Ns), L);
        
        T = Tf*Tb;
%         TWf = TfWf*TbWf; 
        W = (Wf*Wb)';
        
        Rs = P/Ns*eye(Ns);% power allocation embedded in precoding
        Rn = W*W';
        WHT = W*H*T;
%         WHTWf = W*H*TWf; 
        rateTmpEq = log2(det(eye(Ns) + inv(Rn)*WHT*Rs*WHT'));
        rateSpaEq(isnr) = rateSpaEq(isnr) + rateTmpEq;
%         rateTmpWf = log2(det(eye(Ns) + inv(Rn)*WHTWf*Rs*WHTWf'));
%         rateSpaWf(isnr) = rateSpaWf(isnr) + rateTmpWf;
        
        %==================================================================
        % ============ Hybrid RF-baseband precoding =======================
        %==================================================================
%         Tf = zeros(Nt, L);% TX RF precoding
%         Wf = zeros(Nr, L);% RX RF precoding
%         for ichain = 1 : L
%             Tf(:, ichain) = 1/sqrt(Nt)*exp(j*phase(V(:, ichain)));
%             Wf(:, ichain) = 1/sqrt(Nr)*exp(j*phase(U(:, ichain)));
%         end
        
        % =========== PINV ===========
%         Wb = pinv(Wf)*V(:,1:L);
%         Tb = U(:,1:L)'*pinv(Tf);
        % ============================
        
        % ========== SVD =============
        [Ueq, Seq, Veq] = svd(Wf'*H*Tf);
        Tb = Veq(:, 1:Ns);% TX baseband precoding
        Wb = Ueq(:, 1:Ns);% RX baseband precoding
        
%         gain_eq = diag(Seq);% vectorize
%         Gamma_eq = waterfill(P, gain_eq(1:Ns));% power allocation
%         Gamma_eq = diag(Gamma_eq);
%         TbWf = Tb*sqrt(Gamma_eq);

        T = Tf*Tb*sqrt(Ns)/norm(Tf*Tb, 'fro'); 
%         TWf = Tf*TbWf*sqrt(Ns)/norm(Tf*TbWf, 'fro');
        W = (Wf*Wb)';
        
        Rn = W*W'; 
        Rs = P/Ns*eye(Ns);
        WHT = W*H*T;
%         WHTWf = W*H*TWf; 
        rateTmpEq = log2(det(eye(Ns) + inv(Rn)*WHT*Rs*WHT'));
        rateHybEq(isnr) = rateHybEq(isnr) + rateTmpEq;
%         rateTmpWf = log2(det(eye(Ns) + inv(Rn)*WHTWf*Rs*WHTWf'));
%         rateHybWf(isnr) = rateHybWf(isnr) + rateTmpWf;
    end
    isnr
end

rateSvdEq = rateSvdEq/channNum; 
rateHybEq = rateHybEq/channNum; 
rateSpaEq = rateSpaEq/channNum; 
rateSvdWf = rateSvdWf/channNum;
rateHybWf = rateHybWf/channNum;
rateSpaWf = rateSpaWf/channNum;

%%%%%%% Figure plotting for mmWave channel
figure
lw = 1.5;
ms = 6;
plot(SNR, abs(rateSvdEq), 'k-*', 'LineWidth', lw, 'MarkerSize', ms)
hold on
plot(SNR, abs(rateSpaEq), 'b-o', 'LineWidth', lw, 'MarkerSize', ms)
hold on
plot(SNR, abs(rateHybEq), 'r-^', 'LineWidth', lw, 'MarkerSize', ms)
% hold on
% plot(SNR, abs(rateSvdWf), 'k-x', 'LineWidth', lw, 'MarkerSize', ms)
% hold on
% plot(SNR, abs(rateSpaWf), 'b-d', 'LineWidth', lw, 'MarkerSize', ms)
% hold on
% plot(SNR, abs(rateHybWf), 'r-v', 'LineWidth', lw, 'MarkerSize', ms)
hold off
legend('Optimum SVD', 'Sparse precoding & Combining', 'Hybrid', ...
    'Optimum SVD Waterfilling', 'Sparse precoding & Combining Wf', 'Hybrid Wf')
xlabel('SNR (dB)')
ylabel('Sum spectral efficiency (bps/Hz)')
title(sprintf('Nt = %d, Nr = %d, Ncls = %d, Nray = %d, Ns = %d, L = %d',...
    Nt,Nr,Ncls,Nray,Ns,L))
grid
% saveas(gcf, sprintf('mmWaveP2P-Nt%d-Nr%d-Ncls%d-Nray%d-Ns%d-L%d',Nt,Nr,Ncls,...
% %   Nray,Ns,L));% ULA
% saveas(gcf, sprintf('UPAmmWaveP2P-Nt%d-Nr%d-Ncls%d-Nray%d-Ns%d-L%d',Nt,Nr,Ncls,...
%    Nray,Ns,L));% UPA

toc
        
    


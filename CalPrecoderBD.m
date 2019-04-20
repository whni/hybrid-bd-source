%  BD precoder design subfunction r = BD(G, K) --%
%  Input:  H, K*Nr x Nt,  aggregate channel
%          K, User number
%          Ns, # streams per UE
%  Output: r, Nt x K*Ns BD precoder
%          rCombiner, K*Ns x Nr, receiver combiner
%  By Le Liang, UVic, Oct. 27, 2013

function [Precoder, rCombiner] = CalPrecoderBD(H, K, Ns)
if (nargin < 3)
    Ns = size(H,1)/K;
end
[NR, Nt] = size(H);
Nr = NR/K;
F = zeros(Nt, K*Ns);
rCombiner = zeros(K*Ns, Nr);
%%%%%%% Standard BD design based on paper "ZF methods for DL spatial multiplexing..." %%%%%%%%
for ik = 1 : K
    Htmp = H;
    Hk = Htmp(((ik-1)*Nr+1):(ik*Nr),:);
    Htmp(((ik-1)*Nr+1):(ik*Nr),:) = []; % Delete i-user's channel, G_tilde
    [~, ~, V] = svd(Htmp);
    rv = rank(Htmp);
    
    Htmp = Hk*V(:, (rv + 1):Nt);
%     Lr = rank(Htmp);% max #streams per UE
%     Lr = min(Lr, L);
    [Ut, St, Vt] = svd(Htmp);
    F(:, (ik-1)*Ns+1:ik*Ns) = V(:, (rv + 1):Nt)*Vt(: , 1:Ns);
    rCombiner((ik-1)*Ns+1:ik*Ns, :) = Ut(:, 1:Ns)';
end
Precoder = F * sqrt(Ns*K)/ norm(F, 'fro'); % normalize

end
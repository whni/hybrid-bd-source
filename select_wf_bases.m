function [ Wf ] = select_wf_bases( H, K, Ar, Lr )
%SELECT_WF_BASES Summary of this function goes here
%   Detailed explanation goes here
%   Input;
%       H: channel matrix
%       K: number of users
%       Ar: the original ar bases
%       Lr: chain number of UEs
%   Output:
%       Wf: RF combiners

[NR, Nt] = size(H);
Nr = NR/K;
[~, NP] = size(Ar);
Np = NP/K;
Wf = zeros(Nr, Lr*K);
for k = 1:K
    Hk = H((k-1)*Nr+1:k*Nr, :);
    Ark = Ar(:, (k-1)*Np+1:k*Np);
    fnorm = zeros(Np, 1);
    Wff_tmp = zeros(Nr, Lr);
    lr = 0;
    for ip = 1:Np
%         fnorm(ip) = norm(Ark(:, ip)' * Hk, 'fro');
        fnorm(ip) = norm(Ark(:, ip)' * Hk, 1);
    end
    
    th_cor = 1;
    while lr < Lr
        [~, max_pos] = max(fnorm);
        max_cor = max(abs(Ark(:, max_pos)' * Wff_tmp));
        if max_cor <= th_cor
            lr = lr + 1;
            Wff_tmp(:, lr) = Ark(:, max_pos);
%             mat_cor = abs(Wff_tmp' * Wff_tmp);
%             mat_cor(1:(Lr+1):end) = 0;
%             th_cor = max(max(mat_cor)) + 0.1;
%             if th_cor > Lr/10
%                 th_cor = Lr/10;
%             end
        end
%         lr = lr + 1;
%         max_cor = max(abs(Ark' * Wff_tmp), [], 2);
%         [~, max_pos] = max(fnorm .* abs(log2(max_cor)));
%         Wff_tmp(:, lr) = Ark(:, max_pos);
        fnorm(max_pos) = 0;
    end
        
    Wf(:, (k-1)*Lr+1:k*Lr) = Wff_tmp;
end
    
end


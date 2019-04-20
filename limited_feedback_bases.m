function [ Tf, Wf ] = limited_feedback_bases( H, K, At, Ar )
%LIMITED_FEEDBACK_BASES Summary of this function goes here
%   Detailed explanation goes here

[NR, Nt] = size(H);
Nr = NR/K;
[~, NP] = size(Ar);
Np = NP/K;
Wf = zeros(Nr, K);
Tf = zeros(Nt, K);

for k = 1:K
    Hk = H((k-1)*Nr+1:k*Nr, :);
    Ark = Ar(:, (k-1)*Np+1:k*Np);
    Atk = At(:, (k-1)*Np+1:k*Np);
    
    gain_mt = abs(Ark' * Hk * Atk);
    
    max_gain = max(max(gain_mt));
    
    [max_row, max_col] = find(gain_mt == max_gain);
    
    Tf(:, k) = Atk(:, max_col);
    Wf(:, k) = Ark(:, max_row);
end

end


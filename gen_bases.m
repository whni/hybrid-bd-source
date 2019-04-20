function [ at_angle_set ] = gen_bases( H, K, Nb )
%GEN_BASES Summary of this function goes here
%   Detailed explanation goes here
%   Input;
%       H: channel matrix
%       K: number of users
%       Nb: number of bases for each user
%   Output:
%       at_angle_set: the angles of the bases for OMP

[NR, Nt] = size(H);
Nr = NR/K;
d = 0.5;
at_angle_set = zeros(K*Nb, 1);
j = sqrt(-1);
for k = 1:K
    Htmp = H;
    Hk = Htmp((k-1)*Nr+1:k*Nr, :);
    Htmp((k-1)*Nr+1:k*Nr, :) = [];

    [~, S0, V0] = svd(Hk);
    [~, S1, V1] = svd(Htmp);
    L0 = rank(Hk);
    L1 = rank(Htmp);
    num_path = 360;
    at_cor = zeros(num_path, 1);
    for angle = 1:num_path;
        at = exp(j*2*pi*d*sin(angle*2*pi/num_path) *([0:(Nt-1)]') )/sqrt(Nt);
        sum_span0 = sum(abs(V0(:, 1:L0)' * at) .* (diag(S0(1:L0, 1:L0)).^2));
        sum_null0 = sum(abs(V0(:, (L0+1):Nt)' * at));
        
        sum_span1 = sum(abs(V1(:, 1:L1)' * at) .* (diag(S1(1:L1, 1:L1)).^2));
        sum_null1 = sum(abs(V1(:, (L1+1):Nt)' * at));
        at_cor(angle) = (sum_span0 + sum_null1) / (sum_null0 + sum_span1);
    end
    [~, max_pos] = sort(at_cor, 'descend');
    at_angle_set((k-1)*Nb+1:k*Nb) = max_pos(1:Nb)*2*pi/num_path;

end
    
end


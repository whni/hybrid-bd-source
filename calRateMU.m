% r = calRateMU(H, W, T, Rs, NsUE)
% calRateMU calculates achievable rates of DL multi-antenna MU-MIMO systems with TX/RX processing
% Input: 
%   H, NR x Nt, channel instantiation, NR = K x Nr
%   T, Nt x Ns, TX precoder
%   W, Ns x Nr, RX combiner
%   Rs, diagonal Ns x Ns, signal covariance matrix
%   NsUE, 1 x K, stream allocation among users, sum(NsUE) = Ns
% Output: 
%   r, achievable rate of the MU-MIMO system
% By Le Liang, UVic, July 16, 2014


function r = calRateMU(H, T, W, Rs, NsUE)

if (nargin < 5)% default setting, 1 stream & 1 antenna, per user
    NsUE = ones(size(T, 2), 1);
end

Rs = diag(Rs); % vectorize
K = length(NsUE);% #users
Nr = size(H, 1)/K;% #antennas per user, assuming equal
Ns = sum(NsUE);% total #streams

r = 0;% result initialized to zeros

for ik = 1 : K        
    st = sum(NsUE(1:(ik-1)))+1;% start of ik user stream
    ed = sum(NsUE(1:ik));
    
    % ============= Signal power =====================
    Hj = H((Nr*(ik-1)+1):(Nr*ik), :);% user channel of ik
    Tj = T(:, st:ed);% precoder for user ik
    Wj = W(st:ed, :);% combiner for user ik
    Rj = diag(Rs(st:ed));% covariance matrix
    Pj = (Wj*Hj*Tj)*Rj*(Wj*Hj*Tj)'; % signal power term
    
    % ============ interference + noise ==============
    Tintf = T;% precoder for interferers
    Tintf(:, st:ed) = [];
    R_intf = Rs;% covariance of interfering signal
    R_intf(st:ed) = [];
    R_intf = diag(R_intf);
    Pintf = Wj*(eye(Nr)+(Hj*Tintf)*R_intf*(Hj*Tintf)')*Wj';% inter-user interference power term
   
    % ============ rate calculation ==================
    r = r + log2(det(eye(NsUE(ik)) + pinv(Pintf)*Pj));
end




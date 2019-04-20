% Function to calculate achievable rates, optimal ML detectio at the receiver side
% r = calRateOptUE(Incov, H, W, NsUE);
% Input: Incov, input covariance, i.e., power allocation, Ns x Ns
%        H, channel instantialization, (KNr) x Nt
%        W, precoding matrix, column normalized, Nt x Ns
%        NsUE, a vector showing # streams of each UE, K x 1, sum(NsUE) = Ns
% Output: r, achievable rate
% By Le Liang, UVic, Oct. 28, 2013

function r = calRateOptUE(H, W, Incov, NsUE)
if nargin <= 3
    NsUE = ones(size(W,2), 1);% set default, single stream per user
end

Incov = diag(Incov); % vectorize
K = length(NsUE);% # of users
Nr = size(H, 1)/K;% # antenna per user, assuming equal
Ns = sum(NsUE);% total # data streams

rate = 0;

for iK = 1 : K        
    st = sum(NsUE(1:(iK-1)))+1;
    ed = sum(NsUE(1:iK));
    
    %% Signal power term
    Hj = H( (Nr*(iK-1)+1) : (Nr*iK), :);% user channel of i
    Wj = W(:, st:ed);% precoder for user i
    covj = diag(Incov(st:ed));
    Pj = Hj*Wj*covj*(Hj*Wj)'; % signal power term
    
    %% Inter-user interference power term
    Wintf = W;% precoder for interferers
    Wintf(:, st:ed) = [];
    cov_intf = Incov;
    cov_intf(st:ed) = [];
    cov_intf = diag(cov_intf);
    Pintf = (Hj*Wintf)*cov_intf*(Hj*Wintf)';% inter-user interference power term
   
    %% rate Calculation
    rate = rate + log2(det( eye(Nr) + inv( eye(Nr) + Pintf ) * Pj ));
end

r = rate;

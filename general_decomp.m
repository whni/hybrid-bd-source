function [ Frf, Fbb ] = general_decomp( F, Ns, Nrf )
%GENERAL_DECOMP Summary of this function goes here
%   Detailed explanation goes here
%   Input parameter =========
%   F: the original matrix to be decomposed
%   Nrf: the number of chains -- K
%   Output parameter ========
%   Frf: matrix in RF domain
%   Fbb: matrix in baseband

%% Initialization
% deal with zero columns in F (for water-filling)
col_zero = sum(sum(F) == 0);
if col_zero > 0
    a = 1;
end
F = F(:, 1:(Ns - col_zero));
[Nant, ~] = size(F);   

% Use the phase of the first Nrf columns of the objective matrix
j = sqrt(-1);
[U, S, V] = svd(F);
Frf = U*S;
Frf = Frf./abs(Frf);
Frf = [Frf, exp(j*rand(Nant, Nrf - (Ns - col_zero)))];
%Frf = exp(j*rand(Nant, Nrf));
Fbb = ((Frf' * Frf) \ Frf') * F;

% phase_record = zeros(100, 1);
% error_record = zeros(100, 1);
% k = 1;
% phase_record(k) = angle(Frf(1, 5));
% error_record(k) = norm(F - Frf*Fbb, 'fro') / norm(F, 'fro');

%% Iteration
delta_bound = 0.1;        % the upper bound for the increment of the phase in each iteration
pre_err = +inf;         % record the error ||F - F_{RF}F_{BB}||
now_err = norm(F - Frf*Fbb, 'fro') / norm(F, 'fro');
pre_norm = 0;           % record the increment of ||F_{RF}||
now_norm = norm(angle(Frf), 'fro') / sqrt(Nant*Nrf*2*pi);

eps_err = 10e-5;        % the tolerance for the error increment
%tic
while abs(pre_err - now_err) > eps_err
%     k = k + 1;
    % minimize ||F - (Frf_(k)+j{e^(j\phi_{pq})*\delta_{pq}})*Fbb||_F
    Fk = F - Frf*Fbb;
%     Frf_delta = zeros(Nant, Nrf);
%     for p = 1:Nant
%         f_p = Fk(p, :);
%         C_p = j*diag(Frf(p, :))*Fbb;
%         
%         cvx_begin quiet
%             variable d(1, Nrf);
%             minimize(norm(f_p - d*C_p, 'fro'));
%             subject to
%                 abs(d) <= 0.1*ones(1, Nrf);
%         cvx_end
%         Frf_delta(p, :) = d;
%     end
    
    % obtain the next F_{RF}
    cvx_begin quiet
        variable d(Nant, Nrf);
        minimize(norm(Fk - j*(d.*Frf)*Fbb, 'fro'));
        subject to
            abs(d) <= delta_bound*ones(Nant, Nrf);
    cvx_end
    Frf_delta = d;
    Frf = Frf .* exp(j*Frf_delta);
    Fbb = ((Frf' * Frf) \ Frf') * F;
    
    pre_err = now_err;
    now_err = norm(F - Frf*Fbb, 'fro') / norm(F, 'fro');
    
%     phase_record(k) = angle(Frf(1, 5));
%     error_record(k) = norm(F - Frf*Fbb, 'fro') / norm(F, 'fro');
    
    % rescale the bound based on the current error increment
    tang = pre_err - now_err;
    if tang < eps_err*100
        delta_bound = delta_bound * 0.8;
        if delta_bound < 0.1
            delta_bound = 0.1;
        end
    else
        delta_bound = delta_bound * 1.25;
        if delta_bound > 0.5
            delta_bound = 0.5;
        end
    end
    
%     disp(['err = ', num2str(now_err), '; delta_bound = ', num2str(delta_bound)]);
end

disp(['  err = ', num2str(now_err)]);
 
Frf = Frf / sqrt(Nant);
Fbb = [Fbb, zeros(Nrf, col_zero)] * sqrt(Nant);
Fbb = sqrt(Ns) * Fbb / norm(Frf*Fbb, 'fro');
%toc
end


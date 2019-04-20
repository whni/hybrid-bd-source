% Spatially sparse precoding via OMP
% [F_rf, F_bb] = CalSparsePrecoder(A, F_opt, Nch)
% Input: A, candidates for F_rf
%        F_opt, optimal precoder without constant modulus constraints, normzlied
%        frobenius norm
%        Nch, chain number at the BS
% Output: F_bb, baseband precoder
%         F_rf, RF analog preocessing by phase shifters
%         Notes: norm(F_rf*F_bb,'fro')=1;
% By Le Liang, UVic, Oct. 31, 2013

function [F_rf, F_bb] = CalSparsePrecoder(A, F_opt, Nch)

F_rf = [];
F_res = F_opt;

for i = 1 : Nch
    psi = A'*F_res;
    [maxVal, maxInd] = max(diag(psi*psi')); % find the max row magnitude of psi
    F_rf = [F_rf A(:, maxInd)];
    F_bb = pinv(F_rf) * F_opt;
    F_res = F_opt - F_rf*F_bb;
	F_res = F_res/norm(F_res, 'fro');
end
F_bb = sqrt(size(F_opt, 2)) *F_bb/norm(F_rf*F_bb, 'fro');
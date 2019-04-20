% [Wf, Wb] = calMmseCombiner(A, Wmmse, Eyy, L)
% calMmseCombiner computes sparse combiner to approximate Wopt with L chains
% Input:
%   A, Nr x Np, candidates for Wf columns
%   Wmmse, Nr x Ns, uncontrained MMSE combiner
%   Eyy, Nr x Nr, received signal covariance matrix
%   L, scalar, #chains at RX
% Output:
%   Wf, Nr x L, RF combiner
%   Wb, L x Ns, baseband combiner, NOTE: W = (Wf*Wb)';
% By Le Liang, UVic, July 21, 2014

function [Wf, Wb] = calMmseCombiner(A, Wmmse, Eyy, L)
Wf = [];
Wres = Wmmse;
for iL = 1 : L
    phi = A'*Eyy*Wres;
    [maxVal, maxInd] = max(diag(phi*phi'));
    Wf = [Wf A(:, maxInd)];
    Wb = inv(Wf'*Eyy*Wf)*Wf'*Eyy*Wmmse;
    Wres = Wmmse - Wf*Wb;
    Wres = Wres / norm(Wres, 'fro');
end
% --- Generate a clustered channel model for mmWave channels based on "Spatially sparse precoding in millimeter wave MIMO systems"
% --- [H, GainAlpha, Ar, At] = GenChannel (Nt, Nr, Ncls, Nray, sigGain, sigAngle, DirTX, DirRX)
%  Input: 
%      Nt, number of transmit antennas, 1 x 1 for ULA, 1 x 2 for UPA
%      Nr, number of receive antennas, 1 x 1 for ULA, 1 x 2 for UPA
%      Ncls, number of clusters
%      Nray, number of rays per cluster
%      d, normalized antenna spacing
%      sigGain, Ncls x 1, standard deviation of the gain of rays from each cluster, norm(sigGain, 2)^2 = Nt*Nr/Nray, normalized such that E[|H|^2] = Nt*Nr;
%      sigAngle, 4 x 1, if type == 1, reduce to 2 x 1, angular spread of departure azimuth, departure elevation, arrival azimuth, and arrival elevation
%      DirTX, 2 x 2, TX directional antenna parameters    [phi_min, phi_max         phi: azimuth angle, theta: elevation angle
%             if type == 1, reduce to 1 x 2;        theta_min, theta_max]
%      DirRX, 2 x 2, RX directional antenna parameters    [phi_min, phi_max
%             if type == 1, reduce to 1 x 2;        theta_min, theta_max]
%  Output: H, Nr x Nt, composite channel
%      GainAlpha, Ncls*Nray x 1, complex gains of each ray from each cluster
%      Ar, Nr x (Ncls*Nray), aggregate receive array response vector
%      At, Nt x (Ncls*Nray), aggregate transmit array response vector
% example:  Nt = 4, Nr =4, Ncls = 6, Nray = 10
%      sigGain = sqrt(Nt*Nr/(Ncls*Nray))*ones(Ncls, 1); 
%      sigAngle = 7.5*pi/180*ones(4, 1); 
%      DirTX = [-60, 60]*pi/180; DirRX = [-pi, pi];
%      [H, GainAlpha, Ar, At] = GenChannel(Nt, Nr, Ncls, Nray, d, sigGain, sigAngle, DirTX, DirRX)
%      [H, GainAlpha, Ar, At] = GenChannel(Nt, Nr, Ncls, Nray);

% ========  By Le Liang, UVic, Nov. 13, 2013
%           Modified July 7, 2014


function [H, GainAlpha, At, Ar] = GenChannel(Nt, Nr, Ncls, Nray, DirTX, DirRX, d, sigGain, sigAngle)

if (length(Nt) == 1) 
    type = 1; % type 1 for ULA, 2 for UPA
else
    type = 2;
end

% ============== Default Parameters setting ==========================
% ====================================================================
if (nargin  <= 6) % set default
    d = 0.5; % normalized antenna spacing
    if (type == 1) % ULA
        sigGain = sqrt(Nt*Nr/(Ncls*Nray))*ones(Ncls, 1);
        sigAngle = 7.5/180*pi*ones(2, 1);% angular spread
%         DirTX = [-60 60]*pi/180;% omni-directional TX antenna [-pi, pi]
%         DirRX = [-180 180]*pi/180;% omni-directional RX antenna, [-pi, pi]   
    else  % UPA
        sigGain = sqrt(prod(Nt)*prod(Nr)/(Ncls*Nray))*ones(Ncls, 1); % normalized such that E[|H|^2] = Nt*Nr
        sigAngle = 7.5/180*pi*ones(4, 1);
%         DirTX = [-60 60; -60 60]*pi/180; % omnidirectional TX antenna
%         DirRX = [-180 180; -180 180]*pi/180;% omnidirectional RX antenna  
    end
end
% ====================================================================
% ====================================================================


Np = Ncls*Nray; % total #paths
GainAlpha = zeros(Np, 1);% complex gains of paths
j = sqrt(-1);
if type == 1 % ULA antenna
    meanAngle(:, 1) = (2*rand(Ncls, 1)-1)*(DirTX(2)-DirTX(1))/2 + (DirTX(2)+DirTX(1))/2; % mean angle uniformly distributed ~ [-pi, pi]
    meanAngle(:, 2) = (2*rand(Ncls, 1)-1)*(DirRX(2)-DirRX(1))/2 + (DirRX(2)+DirRX(1))/2;
    % col1: mean depart azimuth of each cluster
    % col2: mean arrival azimuth of each cluster
    At = zeros(Nt, Np);
    Ar = zeros(Nr, Np);
    iN = 1;% dummy index for #paths
    while (iN <= Np)
        icls = ceil(iN/Nray);
        sig = sigGain(icls);
        GainAlpha(iN) = sig*(randn(1) + j*randn(1))/sqrt(2);
        
        meanDepart = meanAngle(icls, 1);
        meanArrival = meanAngle(icls, 2);
        sigDepart = sigAngle(1);
        sigArrival = sigAngle(2);
        
        AngD = genLaplacian(meanDepart, sigDepart);
        AngA = genLaplacian(meanArrival, sigArrival);
%         AngD = randLaplaceTrunc(1, 1, meanDepart, sigDepart);% Truncated Laplacian distribution
%         AngA = randLaplaceTrunc(1, 1, meanArrival, sigArrival);
        
        DirGain = (AngD>DirTX(1))*(AngD<DirTX(2))*(AngA>DirRX(1))*(AngA<DirRX(2));
        DirGain = 1;
        if (DirGain == 1)
            Attmp = exp(j*2*pi*d*sin(AngD) *(0:(Nt-1)) )/sqrt(Nt);
            At(:, iN) = Attmp.';
            Artmp = exp(j*2*pi*d*sin(AngA) *(0:(Nr-1)) )/sqrt(Nr);
            Ar(:, iN) = Artmp.';
            iN = iN + 1;        
        end      
    end
    H = Ar * diag(GainAlpha) * At';
 
else % UPA antenna
    % meanAngle = (2*rand(Ncls, 4)-1)*pi; % uniformly distributed mean angles
    meanAngle(:, 1) = (2*rand(Ncls, 1)-1)*(DirTX(1, 2)-DirTX(1, 1))/2 + (DirTX(1, 2)+DirTX(1, 1))/2;% mean departure azimuth of each cluster
    meanAngle(:, 2) = (2*rand(Ncls, 1)-1)*(DirTX(2, 2)-DirTX(2, 1))/2 + (DirTX(2, 2)+DirTX(2, 1))/2;% mean departure elevation of each cluster
    meanAngle(:, 3) = (2*rand(Ncls, 1)-1)*(DirRX(1, 2)-DirRX(1, 1))/2 + (DirRX(1, 2)+DirRX(1, 1))/2;% mean arrival azimuth of each cluster
    meanAngle(:, 4) = (2*rand(Ncls, 1)-1)*(DirRX(2, 2)-DirRX(2, 1))/2 + (DirRX(2, 2)+DirRX(2, 1))/2;% mean arrival elevation of each cluster
    
    Ar = zeros(prod(Nr), Np);
    At = zeros(prod(Nt), Np);
    iN = 1;
    while(iN <= Np)
        icls = ceil(iN/Nray);
        sig = sigGain(icls);
        GainAlpha(iN) = sig * (randn(1) + j*randn(1))/sqrt(2);
        
        meanDepartAzi= meanAngle(icls, 1);
        meanDepartElev = meanAngle(icls, 2);
        meanArrivalAzi = meanAngle(icls, 3);
        meanArrivalElev = meanAngle(icls, 4);
        
        phi_t = genLaplacian(meanDepartAzi, sigAngle(1));        
        theta_t = genLaplacian(meanDepartElev, sigAngle(2));
        phi_r = genLaplacian(meanArrivalAzi, sigAngle(3));
        theta_r = genLaplacian(meanArrivalElev, sigAngle(4));
%         phi_t = randLaplaceTrunc(1, 1, meanDepartAzi, sigAngle(1));
%         theta_t = randLaplaceTrunc(1, 1, meanDepartElev, sigAngle(2));
%         phi_r = randLaplaceTrunc(1, 1, meanArrivalAzi, sigAngle(3));
%         theta_r = randLaplaceTrunc(1, 1, meanArrivalElev, sigAngle(4));
        
        Wt = Nt(1);% TX array width, array dimension: Ht x Wt
        Ht = Nt(2);% TX array height
        Wr = Nr(1);% RX array width
        Hr = Nr(2);% RX array height
        
        DirGain = (phi_t > DirTX(1,1))*( phi_t < DirTX(1,2))*(theta_t > DirTX(2,1))*(theta_t < DirTX(2,2))* ...
        (phi_r > DirRX(1,1))*(phi_r < DirRX(1,2))*(theta_r > DirRX(2,1))*(theta_r < DirRX (2,2));   
        DirGain = 1;
        if (DirGain == 1)
            for it = 1 : Wt*Ht
                iH = floor((it-1)/Wt);
                iW = it - iH*Wt - 1;
                At(it, iN) = 1/sqrt(Wt*Ht) * exp(j*2*pi*d*( iW*sin(phi_t)*sin(theta_t) + iH*cos(theta_t) ));
            end

            for ir = 1 : Wr*Hr
                iH = floor((ir-1)/Wr);
                iW = ir - iH*Wr - 1;
                Ar (ir, iN) = 1/sqrt(Wr*Hr) * exp(j*2*pi*d*( iW*sin(phi_r)*sin(theta_r) + iH*cos(theta_r) ));
            end
            iN = iN + 1;
        end
    end
    H = Ar * diag(GainAlpha) * At';
end
        

%%%%%%%%%% Generate a Laplacian variable with parameter {mu, sig}
%   mean = mu, std deviation = sig
%   Reference: wikipedia: http://en.wikipedia.org/wiki/Laplace_distribution
function x = genLaplacian(mu, sig) 
b = sig/sqrt(2);
u = rand - 1/2;
x = mu - b*sign(u)*log(1-2*abs(u))/log(exp(1));


%RANDLAPLACETRUNC generate i.i.d. laplacian random number drawn from
%truncated laplacian distribution over interval [-pi,pi) with mean mu and
%standard deviation sigma.
%   Parameters:
%   mu      : mean
%   sigma   : standard deviation
%   [m, n]  : the dimension of y.
%   Default mu = 0, sigma = 1.
%
%   Reference:
%   [1] http://en.wikipedia.org./wiki/Laplace_distribution.
%   [2] Antonio Forenza, David J. Love, and Robert W. Heath Jr.,
%   "Simplified spatial correlation models for clustered MIMO channels with
%   different array configurations", IEEE Transaction on Vehicular Technology,
%   vol. 56, no. 4, July 2007.
%
%   Copyright 2014 Leyuan Pan, Xiaodai Dong
%   Designed by Leyuan Pan (leyuanpan@gmail.com)
%   $Revision: 0.1 $  $Date: 2014/06/24 14:05 $

function y  = randLaplaceTrunc(m, n, mu, sigma)
b = sigma/sqrt(2);
u = 2*(1-exp(-pi/b))*(rand(m, n)-0.5);
y = mu-b*sign(u).*log(1-abs(u));




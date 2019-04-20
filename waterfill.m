% waterfilling power allocation
% Input: P, total power
%        h, n x 1, channel gain/noise
% Output: r, n x 1, power allocation result
%		  s, water surface
% By Le Liang, UVic, Feb. 21, 2014

function [r s]= waterfill(P, h)

bottom = 1./(abs(h)).^2; % bottom of vessel
n = length(h);
bsort = sort(bottom);% sort bottom in ascending order

m = 1;
if (n > 1)
    pt = m*bsort(m+1) - sum(bsort(1:m));% water level staircase
    while (pt < P && m < (n-1))
        m = m + 1;
        pt = m*bsort(m+1) - sum(bsort(1:m));
    end
    if (m == (n-1))% water level higher than the top
        m = n;
    end
end
surface = (P + sum(bsort(1:m)))/m; % water surface
r = max((surface - bottom), 0);

s = surface;
end


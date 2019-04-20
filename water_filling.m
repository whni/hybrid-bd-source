function [ vec_power_alloc ] = water_filling( signal_power, vec_base_level )
%WATER_FILLING Summary of this function goes here
%   Input parameters =================
%   signal_power: the amount of power that is available for the transmitter
%   vec_base_level: the vector that constains the base levels for the water-filling
%   Output parameters ================
%   vec_power_alloc: the vector that shows the water-filling power allocation

%   grab the basic information of the data streams
num_stream = length(vec_base_level);
num_alloc = num_stream;
tmp_base_level = vec_base_level;

%   search the optimal water level
water_level = (sum(tmp_base_level) + signal_power) / num_alloc;
[max_level, max_pos] = max(tmp_base_level);
while water_level < max_level
    tmp_base_level(max_pos) = [];
    num_alloc = num_alloc - 1;
    water_level = (sum(tmp_base_level) + signal_power) / num_alloc;
    [max_level, max_pos] = max(tmp_base_level);
end

%   do water filling based on the water level
vec_power_alloc = zeros(num_stream, 1);
for cnt = 1:num_stream
    vec_power_alloc(cnt) = max(water_level - vec_base_level(cnt), 0);
end

end


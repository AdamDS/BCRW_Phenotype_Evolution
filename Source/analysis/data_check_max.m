%% data_check_max.m
%[cout] = data_check_max(data,compval)
%

function [cout] = data_check_max(data,compval), 
cout = zeros(2,1);
cout(1,1) = max(max(data));
if cout(1,1)<=compval, cout(2,1) = 1; end
end
%% data_check_min.m
%[cout] = data_check_min(data,compval)
%

function [cout] = data_check_min(data,compval), 
cout = zeros(2,1);
cout(1,1) = min(min(data));
if cout(1,1)>=compval, cout(2,1) = 1; end
end
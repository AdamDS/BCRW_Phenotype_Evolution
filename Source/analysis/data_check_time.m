%% data_check_time.m
%[cout] = data_check_time(data,compval)
%

function [cout] = data_check_time(data,compval), 
cout = zeros(2,1);
cout(1,1) = length(data);
find_ = find(data>=compval);
if cout(1,1)==length(find_), cout(1,1) = 1; end
end
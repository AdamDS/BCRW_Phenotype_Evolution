%% data_check_length.m
% [cout] = data_check_length(data,compval)
%

function [cout] = data_check_length(data,compval), 
cout = zeros(2,1);
cout(1,1) = length(data);
if cout(1,1)<=compval, cout(2,1) = 1; end
end
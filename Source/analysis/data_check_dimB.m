%% data_check_dimB.m
%[cout] = data_check_dimB(data,compval)
%

function [cout] = data_check_dimB(data,compval), 
cout = zeros(2,1);
cout(1,1) = size(data,2);
if compval(1)==cout(1,1), cout(2,1) = 1;  end
end
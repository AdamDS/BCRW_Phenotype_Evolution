%% data_check_dimA.m
%[cout] = data_check_dimA(data,compval)
%

function [cout] = data_check_dimA(data,compval), 
cout = zeros(2,1);
cout(1,1) = size(data,1);
if compval(1)==cout(1,1), cout(2,1) = 1;  end
end
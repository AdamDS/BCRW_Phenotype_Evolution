%% data_check_initial.m
%[cout] = data_check_initial(data,compval)
%

function [cout] = data_check_initial(data,compval), 
cout = zeros(2,1);
cout(1,1) = data(1);
if compval==cout(1,1), cout(2,1) = 1;  end
end
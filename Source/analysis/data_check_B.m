%% data_check_B.m
% [cout,str_out] = data_check_B(flag,str_data,data,compval,comparisons)
% Checks length, min, & max of data.

function [cout,str_out] = data_check_B(flag,str_data,data,compval,comparisons), 
global str_min str_max str_length;
cout = zeros(2,3);
str_out = {};
if length(compval)~=3,  
  str_out = {['MISSING CHECK DATA: ' str_data '\n']};
else, 
  if flag==1, 
    i = 1;
    cout(:,i) = data_check_length(data,compval(i));
    str_out(i,:) = print_check_result(str_data,str_length,cout(1,i),comparisons(i),compval(i),cout(2,i));
    i = 2;
    cout(:,i) = data_check_min(data,compval(i));
    str_out(i,:) = print_check_result(str_data,str_min,cout(1,i),comparisons(i),compval(i),cout(2,i));
    i = 3;
    cout(:,i) = data_check_max(data,compval(i));
    str_out(i,:) = print_check_result(str_data,str_max,cout(1,i),comparisons(i),compval(i),cout(2,i));
  else, 
    str_out = {['FAIL: ' str_data ' DNE \n']};
  end
end
end

% function [cout] = data_check_B(flag,data,compval), 
% cout = zeros(2,3);
% if flag==1, 
%   cout(:,1) = data_check_length(data,compval(1));
% 
%   cout(:,2) = data_check_min(data,compval(2));
% 
%   cout(:,3) = data_check_max(data,compval(3));
% end
% end
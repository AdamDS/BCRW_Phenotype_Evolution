%% data_check_A.m
% [cout,str_out] = data_check_A(flag,str_data,data,compval,comparisons)
% Checks initial, min, and length values of data.

function [cout,str_out] = data_check_A(flag,str_data,data,compval,comparisons), 
global str_min str_length str_ini;
cout = zeros(2,3);
str_out = {};
if flag==1, 
  i = 1;
  cout(:,i) = data_check_initial(data,compval(i));
  str_out(i,:) = print_check_result(str_data,str_ini,cout(1,i),comparisons(i),compval(i),cout(2,i));
  i = 2;
  cout(:,i) = data_check_min(data,compval(i));
  str_out(i,:) = print_check_result(str_data,str_min,cout(1,i),comparisons(i),compval(i),cout(2,i));
  i = 3;
  cout(:,i) = data_check_length(data,compval(i));
  str_out(i,:) = print_check_result(str_data,str_length,cout(1,i),comparisons(i),compval(i),cout(2,i));
else, 
  str_out = {['FAIL: ' str_data ' DNE \n']};
end
end





% cout(:,1) = data_check_initial(data,compval(1));
% 
% cout(:,2) = data_check_min(data,compval(2));
% 
% cout(:,3) = data_check_time(data,compval(3));
% end
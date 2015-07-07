%% data_check_C.m
% [cout,str_out] = data_check_C(flag,str_data,data,compval,comparisons)
% Checks rows, cols, min, & max of data.

function [cout,str_out] = data_check_C(flag,str_data,data,compval,comparisons), 
global str_min str_max str_rows str_cols;
cout = zeros(2,4);
str_out = {};
if length(compval)~=4,  
  str_out = {['MISSING CHECK DATA: ' str_data '\n']};
else, 
  if flag==1, 
    i = 1;
    cout(:,i) = data_check_dimA(data,compval(i));
    str_out(i,:) = print_check_result(str_data,str_rows,cout(1,i),comparisons(i),compval(i),cout(2,i));
    i = 2;
    cout(:,i) = data_check_dimB(data,compval(i));
    str_out(i,:) = print_check_result(str_data,str_cols,cout(1,i),comparisons(i),compval(i),cout(2,i));
    i = 3;
    cout(:,i) = data_check_min(data,compval(i));
    str_out(i,:) = print_check_result(str_data,str_min,cout(1,i),comparisons(i),compval(i),cout(2,i));
    i = 4;
    cout(:,i) = data_check_max(data,compval(i));
    str_out(i,:) = print_check_result(str_data,str_max,cout(1,i),comparisons(i),compval(i),cout(2,i));
  else, 
    str_out = {['FAIL: ' str_data ' DNE \n']};
  end
end
end


% function [cout] = data_check_C(flag,data,compval), 
% cout = zeros(2,4);
% if flag==1, 
%   cout(:,1) = data_check_dimA(data,compval(1));
%   
%   cout(:,2) = data_check_dimB(data,compval(2));
% 
%   cout(:,2) = data_check_min(data,compval(3));
% 
%   cout(:,3) = data_check_max(data,compval(4));
% end
% end

% row_tcs,col_tcs] = size(trace_cluster_seed); 
% if sum_p==row_tcs,  c_rowtcs = 1; end 
% if (limit-1)==col_tcs,  c_coltcs = 1; end 
% min_tcs = min(min(trace_cluster_seed));
% if 1<=min_tcs,  c_mintcs = 1; end
% max_tcs = max(max(trace_cluster_seed));
% if max_p>=max_tcs,  c_maxtcs = 1; end
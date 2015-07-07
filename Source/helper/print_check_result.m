%% print_check_result.m
% [str_out] = print_check_result(str_value,str_type,value,str_compare,measure,check)
% Generates check result based on the data (str_value) desired value 
% (str_data) compared to (str_compare) the measure of the data (measure). 
%

function [str_out] = print_check_result(str_value,str_type,value,str_compare,measure,check), 
global s;
if exist('measure','var'), 
  if check, 
    str_out = {['pass: ' str_value s str_type s num2str(value) s ...
                       str_compare s num2str(measure) '\n']};
  else, 
    str_out = {['FAIL: ' str_value s str_type s num2str(value) s ...
                       '!' str_compare s num2str(measure) '\n']};
  end
else, 
  str_out = {[]};
end
end
%% exist_load.m
% [data,flag,error] = exist_load(data_name,flag,do_print)
% This function will attempt to load the data file passed to it, data_name.
% If successful, the data is loaded.
% Otherwise, an error message is displayed and the flag, which is boolean,
% will toggle to the opposing state of it's input value.
%
% This is a modified variation of try_catch_load (and even uses that
% function), but it makes use of the function mat_exist. If data_name
% exists, then the data is loaded. If not, then the error is obtained by
% going through try_catch_load. The reason for this change is that
% try_catch_load tends to bog down execution times and can cause matlab to
% seemingly hang on the function. If there is an error, only then will
% there be a problem.
%
% -ADS 5*2*13

function [data,flag,error] = exist_load(data_name,flag,do_print), 
error = []; data = [];
if ~exist('flag','var'), flag = 1;  end
if ~exist('do_print','var'), do_print = 0;  end

if mat_exist(data_name), 
  [data,flag,error] = try_catch_load(data_name,flag,do_print);
else, 
%   data = load(data_name);
  flag = 0;
end
end
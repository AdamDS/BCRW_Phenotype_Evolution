%% try_catch_load.m
% [data,flag,error] = try_catch_load(data_name,flag,do_print)
% This function will attempt to load the data file passed to it, data_name.
% If successful, the data is loaded.
% Otherwise, an error message is displayed and the flag, which is boolean,
% will toggle to the opposing state of it's input value.

function [data,flag,error] = try_catch_load(data_name,flag,do_print), 
error = []; data = [];
if ~exist('flag','var'), flag = 1;  end
if ~exist('do_print','var'), do_print = 0;  end
found0 = 0; found1 = 0; found2 = 0; found3 = 0; tried = 0;
% need to try four variations if they exist
try,  
  data = load(data_name);
  found0 = 1;
catch error,
  found0 = 0;
  tried = tried +1;
end
try,  
  data_name2 = get_alt_split_name(data_name);
  data = load(data_name2);
  found1 = 1;
catch error,  
  found1 = 0;
  tried = tried +1;
end
try,  
  data_name3 = get_alt_sort_name(data_name);
  data = load(data_name3);
  found2 = 1;
catch error,  
  found2 = 0;
  tried = tried +1;
end
try,  
  data_name4 = get_alt_splitsort_name(data_name);
  data = load(data_name4);
  found3 = 1;
catch error,  
  found3 = 0;
  tried = tried +1;
end
if ~found0 && ~found1 && ~found2 && ~found3,  
  if do_print, fprintf('Tried %d %s \n',tried,error.message);  end
  flag = mod(flag+1,2);
end
end %function
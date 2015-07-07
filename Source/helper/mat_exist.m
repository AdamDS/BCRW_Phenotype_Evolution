%% mat_exist.m
% [result] = mat_exist(data_name)
% Boolean check (result) on whether data_name exists.
% There is a check for the extension '.mat'. If this check fails, then 
% the exist check applies '.mat' to the input data_name. 
%
% Care should be taken when inputing data_name both in following naming 
% schemes with spelling and spacings or underscores & in the extension 
% applied.
%
% -ADS 8*27*13
%

function [result] = mat_exist(data_name),  
global SIMOPTS;
% initialize function output
result = 0;
%% Determine if data_name includes the appropriate extension '.mat'
% answers if data_name has the extension '.mat'
has_extension = 0;
% finds the dots
dot = find(data_name=='.');

if dot, %if there is a dot
  %looks for instances of 'm', 'a', and 't'
  ms_as_ts = find(data_name=='m' | data_name=='a' | data_name=='t');
  if ms_as_ts,  %if an 'm', 'a', or 't' exists
    mat = find(ms_as_ts>dot(end));  %get instances that occur after last '.'
    if numel(mat)==3 && ...  %if mat has the right numel
       mat(1)+1==mat(2) && mat(1)+2==mat(3) && ...  %if m is before a and m is two before t
       ms_as_ts(mat(3))==length(data_name), %if mat is at the very end of data_name
      has_extension = 1;
    end
  end
end

if has_extension, 
  %*\Data\*
  ex = exist([data_name],'file');
  %*\Data\#_#\*
  data_name2 = get_alt_split_name(data_name);
  ex2 = exist([data_name2],'file');
  %*\Data\needs_sorting\*
  data_name3 = get_alt_sort_name(data_name);
  ex3 = exist([data_name],'file');
  %*\Data\needs_sorting\#_#\*
  data_name4 = get_alt_splitsort_name(data_name);
  ex4 = exist([data_name4],'file');
else, 
  %try *\Data\*
  ex = exist([data_name '.mat'],'file');
  %*\Data\#_#\*
  data_name2 = get_alt_split_name(data_name);
  ex2 = exist([data_name2 '.mat'],'file');
  %*\Data\needs_sorting\*
  data_name3 = get_alt_sort_name(data_name);
  ex3 = exist([data_name3 '.mat'],'file');
  %*\Data\needs_sorting\#_#\*
  data_name4 = get_alt_splitsort_name(data_name);
  ex4 = exist([data_name4 '.mat'],'file');
end

if ex==2 || ex2==2 || ex3==2 || ex4==2, 
    result = 1; 
else, 
  result = 0; 
end
end


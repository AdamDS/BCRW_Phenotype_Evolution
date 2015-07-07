%% get_alt_split.m
% [alt_name] = get_alt_split_name(data_name)
%
% Adds or removes #_#\ from directory name at the beginning of data_name.
% If data_name has terminal directory as:
% \Data\needs_sorting\#_#\, then alt_name dir ends as \Data\needs_sorting\
% \Data\#_#\, then alt_name dir ends as \Data\
% \Data\needs_sorting\, then alt_name dir ends as \Data\needs_sorting\#_#\
% \Data\, then alt_name directory ends as \Data\#_#
%
% NOTE: It does not matter whether needs_sorting is in the directory part
% of data_name or not because #_# is always at the end of the directory.
% -ADS 8*27*13

function [alt_name] = get_alt_split_name(data_name),  
global SIMOPTS;
alt_name = data_name;
enddir = find(data_name=='\' | data_name=='/'); %get end of dir
if SIMOPTS.split,  %has *\#_#\* in name, so cut it
  if SIMOPTS.for_external, %is *\Data\needs_sorting\#_#\*
    alt_name = data_name([1:enddir(end-1) (enddir(end)+1):end]); %\Data\needs_sorting\
  else, %is *\Data\#_#\*
    alt_name = data_name([1:enddir(end-1) (enddir(end)+1):end]); %\Data\
  end
else, %does not include *\#_#\* in name, so add it
  guess = split_cd(data_name(1:enddir(end)),SIMOPTS.run,10,0,0); %guess split = 10
  if SIMOPTS.for_external, %is *\Data\needs_sorting\*
    alt_name = [guess data_name((enddir(end)+1):end)]; %\Data\needs_sorting\#_#\
  else, %is *\Data\*
    alt_name = [guess data_name((enddir(end)+1):end)]; %\Data\#_#\
  end
end
end
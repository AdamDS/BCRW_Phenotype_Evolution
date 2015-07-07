%% get_alt_sort.m
% [alt_name] = get_alt_sort_name(data_name)
%
% Adds or removes needs_sorting\ from directory name at the beginning of data_name.
% If data_name has terminal directory as:
% \Data\needs_sorting\#_#\, then alt_name dir ends as \Data\#_#\
% \Data\needs_sorting\, then alt_name dir ends as \Data\
% \Data\#_#\, then alt_name dir ends as \Data\needs_sorting\#_#\
% \Data\, then alt_name directory ends as \Data\needs_sorting\
%
% -ADS 8*27*13

function [alt_name] = get_alt_sort_name(data_name),  
global SIMOPTS;
alt_name = data_name;
enddir = find(data_name=='\' | data_name=='/'); %get end of dir
if SIMOPTS.for_external,  %has *\needs_sorting\* in name, so cut it
  if SIMOPTS.split, %is *\Data\needs_sorting\#_#\*
    alt_name = data_name([1:enddir(end-2) (enddir(end-1)+1):end]); %\Data\#_#\
  else, %is *\Data\needs_sorting\*
    alt_name = [data_name(1:enddir(end-1)) data_name((enddir(end)+1):end)]; %\Data\
  end
else, %does not include *\needs_sorting\* in name, so add it
  if SIMOPTS.split, %is *\Data\#_#\*
    alt_name = [data_name(1:enddir(end-1)) 'needs_sorting' ...
      data_name(enddir(end-1):end)]; %\Data\needs_sorting\#_#\
  else, %is *\Data\*
    alt_name = [data_name(1:enddir(end)) 'needs_sorting' ...
      data_name(enddir(end):end)]; %\Data\needs_sorting\
  end
end
end
%% get_alt_splitsort_name.m
% [alt_name] = get_alt_splitsort_name(data_name)
%
% Adds or removes needs_sorting\ & #_#\ from directory name at the beginning of data_name.
% If data_name has terminal directory as:
% \Data\needs_sorting\#_#\, then alt_name dir ends as \Data\
% \Data\needs_sorting\, then alt_name dir ends as \Data\#_#\
% \Data\#_#\, then alt_name dir ends as \Data\needs_sorting\
% \Data\, then alt_name directory ends as \Data\needs_sorting\#_#\
%
% -ADS 8*28*13

function [alt_name] = get_alt_splitsort_name(data_name),  
global SIMOPTS;
alt_name = data_name;
enddir = find(data_name=='\' | data_name=='/'); %get end of dir
if SIMOPTS.for_external,  %has *\needs_sorting\* in name, so cut it
  if SIMOPTS.split, %is *\Data\needs_sorting\#_#\*
    alt_name = data_name([1:enddir(end-2) (enddir(end)+1):end]); %\Data\
  else, %is *\Data\needs_sorting\*
    guess = split_cd(data_name(1:enddir(end-1)),SIMOPTS.run,10,0,0); %guess split = 10
    alt_name = [guess data_name((enddir(end)+1):end)]; %\Data\#_#\
  end
else, %does not include *\needs_sorting\* in name, so add it
  if SIMOPTS.split, %is *\Data\#_#\*
    alt_name = [data_name(1:enddir(end-1)) 'needs_sorting' ...
      data_name(enddir(end):end)]; %\Data\needs_sorting\
  else, %is *\Data\*
    guess = split_cd(data_name(1:enddir(end)),SIMOPTS.run,10,0,0); %guess split = 10
    slashes = find(guess=='\' | guess=='/');
    alt_name = [guess(1:slashes(end-1)) 'needs_sorting' guess(slashes(end-1):end) ...
      data_name((enddir(end)+1):end)]; %\Data\needs_sorting\#_#\
  end
end
end
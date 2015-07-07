function [alt_name] = get_alternative_name(data_name),  
global SIMOPTS;
alt_name = data_name;
enddir = find(data_name=='\' | data_name=='/'); %get end of dir
if SIMOPTS.for_external,  %has *\needs_sorting\* in name, so cut it
  if SIMOPTS.split, %is *\needs_sorting\#_#\*
    alt_name = data_name([1:enddir(end-2) (enddir(end-1)+1):end]); %checked, pass
  else, %is *\needs_sorting\*
    alt_name = data_name([1:enddir(end-1) (enddir(end)+1):end]);
  end
else, %does not include *\needs_sorting\* in name, so add it
  if SIMOPTS.split, %is *\Data\#_#\*
    alt_name = [data_name(1:enddir(end-1)) 'needs_sorting\' ...
      data_name(enddir(end-1):end)];
  else, %is *\Data\*
    alt_name = [data_name(1:enddir(end)) 'needs_sorting\' ...
      data_name(enddir(end):end)];
  end
end
end
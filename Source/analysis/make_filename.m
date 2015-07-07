function [filename] = make_filename(data_name,base_name,run_name);
global SIMOPTS;
filename = [];
if SIMOPTS.limit==3
  filename = [data_name '_' generalize_base_name(base_name)];
else
  filename = [data_name '_' generalize_base_name(base_name) '_' int2str(SIMOPTS.limit) '_limit'];
end
end
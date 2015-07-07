global SIMOPTS;
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
if only_lt==0
  filename = generalize_base_name(base_name);
  lbn = length(mutability); ldm = length(death_max);  sims_length = length(SIMS);
  if lbn>=ldm,  d = lbn;  else, d = ldm;  end
  headers = cell(1,sims_length +3);

  headers(1,1) = cellstr('Mutability');
  headers(1,2) = cellstr('Avg');
  headers(1,3) = cellstr('Std'); 
  for r = 4:(sims_length +3)
    headers(1,r) = cellstr(['Run ' int2str(r-3)]);
  end
  this_alpha = ceil(((sims_length+3)/26)-1);
  if this_alpha==0, end_alpha = [alphabet(mod(sims_length+3,length(alphabet)))];
  else, end_alpha = [alphabet(this_alpha) alphabet(mod(sims_length+3,length(alphabet)))]; end
  part1 = ['A1:' end_alpha '1'];
  part2 = ['A2:' end_alpha int2str(d +size(headers,1))];
else
  get_lifetimes;
  lbn = length(mutability); ldm = length(death_max);  sims_length = length(SIMS);
  filename = make_data_name('lifetimes',base_name,'',1);
  if lbn>=ldm,  d = lbn;  else, d = ldm;  end
  headers = cell(d +1,1);
  sheet_num = 1;
  if lbn<ldm
    headers(1,1) = cellstr('Random Death');
  else
    headers(1,1) = cellstr('Mutability');
  end
  for c = 2:(sims_length +1)
    headers(1,c) = cellstr(['Run ' int2str(c-1)]);
  end
  this_alpha = ceil(((sims_length+1)/26)-1);
  if this_alpha==0, end_alpha = [alphabet(mod(sims_length+1,length(alphabet)))];
  else, end_alpha = [alphabet(this_alpha) alphabet(mod(sims_length+1,length(alphabet)))]; end
  part1 = ['A1:' end_alpha '1'];
  part2 = ['A2:' end_alpha int2str(d +size(headers,1))];
end

here = cd; %where the data is
cd(output); %where you want to save the exported data

sheet_num = ['Lifetimes'];
xlswrite(filename,headers,sheet_num,part1);
if lbn>=ldm
  xlswrite(filename,[mutability' AVG_NGENS STD_NGENS ngens],sheet_num,part2);
else
  xlswrite(filename,[death_max' AVG_NGENS STD_NGENS ngens],sheet_num,part2);
end
cd(here)
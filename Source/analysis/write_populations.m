global SIMOPTS;
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

filename = generalize_base_name(base_name);
lbn = length(mutability); sims_length = length(SIMS);
if by_generation==0, 
  sheet_num = ['Populations'];
  headers = cell(1,sims_length +3);
  headers(1,1) = cellstr('Mutability'); headers(1,2) = cellstr('Avg');  headers(1,3) = cellstr('Std'); 
  for r = 1:sims_length
    headers(1,r+3) = cellstr(['Run ' int2str(r)]);
  end
  this_alpha = ceil(((sims_length+3)/26)-1);
  if this_alpha==0, end_alpha = [alphabet(mod(sims_length+3,length(alphabet)))];
  else, end_alpha = [alphabet(this_alpha) alphabet(mod(sims_length+3,length(alphabet)))]; end
  part1 = ['A1:' end_alpha int2str(size(headers,1))]; %for the headers
  part2 = ['A' int2str(size(headers,1)+1) ':' end_alpha int2str(lbn +size(headers,1))]; %for the data
else, 
  sheet_num = ['Populations v Generations']
  headers = cell(1,1);
  headers(1,1) = cellstr('Mutability');
  this_alpha = ceil(((lbn+1)/26)-1);
  if this_alpha==0, end_alpha = [alphabet(mod(lbn+1,length(alphabet)))];
  else, end_alpha = [alphabet(this_alpha) alphabet(mod(lbn+1,length(alphabet)))]; end
  part1 = ['A1:A' int2str(size(headers,1))];
  part2 = ['B1:' end_alpha int2str(NGEN +size(headers,1))];
end

here = cd;  %what directory this is in now
cd(output); %where you selected to save the data

if by_generation==0
  xlswrite(filename,headers,sheet_num,part1);
  xlswrite(filename,[mutability' AVG_POPS STD_POPS avg_pop],sheet_num,part2);
elseif by_generation==1 && by_window==0
  xlswrite(filename,headers,sheet_num,part1);
  xlswrite(filename,[1:NGEN]',sheet_num,['A3:A' int2str(NGEN+2)]);
  xlswrite(filename,[mutability; AVG_GEN_POPS'],sheet_num,part2);  
elseif by_generation==1 && by_window==1
  
end
cd(here)  %return to the directory this started
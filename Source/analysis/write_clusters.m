global SIMOPTS;
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
filename = generalize_base_name(base_name);
lbn = length(mutability); sims_length = length(SIMS);
headers = cell(1,sims_length +3);

headers(1,1) = cellstr('Mutability'); 
headers(1,2) = cellstr('Avg');
headers(1,3) = cellstr('Std'); 
for r = 4:(sims_length +3)
  headers(1,r) = cellstr(['Run ' int2str(r-3)]);
end

here = cd;
cd(output);

this_alpha = ceil(((sims_length+3)/26)-1);
if this_alpha==0, end_alpha = [alphabet(mod(sims_length+3,length(alphabet)))];
else, end_alpha = [alphabet(this_alpha) alphabet(mod(sims_length+3,length(alphabet)))]; end
part1 = ['A1:' end_alpha int2str(size(headers,1))];
part2 = ['A' int2str(size(headers,1)+1) ':' end_alpha int2str(lbn +size(headers,1))];

sheet_num = ['Clusters'];
xlswrite(filename,headers,sheet_num,part1);
xlswrite(filename,[mutability' AVG_NCS STD_NCS avg_nc],sheet_num,part2);
cd(here)
global SIMOPTS;
alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

filename = [generalize_base_name(base_name)];
lbn = length(mutability); sims_length = length(SIMS);
headers = cell(1,sims_length +3);

%% R
headers(1,1) = cellstr('Mutability');
headers(1,2) = cellstr('Avg');
headers(1,3) = cellstr('Std');
for r = SIMS
  headers(1,r+3) = cellstr(['Run ' int2str(r)]);
end
part1 = ['A1:' alphabet(sims_length +3) int2str(size(headers,1))];
part2 = ['A' int2str(size(headers,1)+1) ':' alphabet(sims_length +3) int2str(size(headers,1) +lbn)];
here = cd;
cd(output);

sheet_num = ['R'];
xlswrite(filename,headers,sheet_num,part1);
xlswrite(filename,[mutability' AVG_R STD_R avgR],sheet_num,part2);

%% c
sheet_num = 'c';
xlswrite(filename,headers,sheet_num,part1);
xlswrite(filename,[mutability' AVG_c STD_c avgc],sheet_num,part2);
cd(here)
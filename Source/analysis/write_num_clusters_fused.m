global SIMOPTS;
sheet_num = 1;
filename = [make_data_name('num_clusters_fused',base_name,'',0)];
lbn = length(mutability); sims_length = length(SIMS);
runn = cell(1,sims_length);
headers = cell(2,sims_length +1);
for r = 1:sims_length
  runn(1,r) = cellstr(['Run ' int2str(r)]);
end

headers(1,1) = cellstr('# fused');   
headers(2,1) = cellstr('Mutability'); headers(2,2:(size(headers,2))) = runn;
headers(2,sims_length +2) = cellstr('Avg');
headers(2,sims_length +3) = cellstr('Std');

here = cd;
cd(output);

alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
xlswrite(filename,headers,sheet_num,['A1:' alphabet(sims_length +3) '2']);
xlswrite(filename,[mutability' avgR AVG_R STD_R],sheet_num,...
  ['A3:' alphabet(sims_length +3) int2str(2 +lbn)]);
cd(here)
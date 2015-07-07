here = cd;
cd 'C:\Users\Amviot\Desktop'
filename = ['Whittaker_' generalize_base_name(base_name,exp_type,reproduction)];
u = 0;  v = 0;
sdt = 0;
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
for i = 1:2:2*length(bn)
  I = ceil(i/2);
  u = v +1; v = sum(num_rel_abundances(1:i));
  this_set = [RELATIVE_ABUNDANCES(:,u:v)'];
  group = floor((i-1)/26);
  this_letter = mod((i-1),26) +1;
  last = int2str(length(this_set));
  if i<=26
    xlswrite(filename,[bn(I); this_set(:,1)],[letters(this_letter) int2str(1) ':' ...
                                letters(this_letter) num2str(last)]);
    xlswrite(filename,[bn(I); this_set(:,2)],[letters(this_letter+1) int2str(1) ':' ...
                                letters(this_letter+1) num2str(last)]);
  else
    xlswrite(filename,[bn(I); this_set(:,1)],[letters(group) letters(this_letter) int2str(1) ':' ...
                                letters(group) letters(this_letter) num2str(last)]);
    xlswrite(filename,[bn(I); this_set(:,2)],[letters(group) letters(this_letter+1) int2str(1) ':' ...
                                letters(group) letters(this_letter+1) num2str(last)]);
  end
end
cd(here)
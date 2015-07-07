% function [data_dir] = get_Data_dir(data_name)
data_dir = [];
data_name = 'C:\Babies_Root\Bacterial\Mu_x\Uniform\2_Flatscape\2x2_basic_map_size\3_IPOP\2_NGEN\Data\';
D = find(data_name=='D');
a = find(data_name=='a');
t = find(data_name=='t');
if length(D) && length(a) && length(t), 
  for i = D,  
    for j = a,  
      for k = t,  
        [[i j k];[i j-1 k-2]]
        if i==(j-1) && i==(k-2),  
          data_dir = data_name(1:k+2);
        end
      end
    end
  end
end
% end
%% print_report
% [] = print_report(outs)
% prints each line of output results
%
function [] = print_report(outs), 
I = size(outs,1);
if I==1, 
  fprintf(outs{1});
else, 
  for i = 1:I,  
    J = size(outs{i},2);  
    for j = 1:J,  
      %if outs{i}{j}=='F',   
        fprintf(outs{i}{j});
      %end
    end
  end
end
end
%% print_report_pass
% [] = print_report_pass(outs)
% prints each line of output results
%

function [] = print_report_pass(outs), 
I = size(outs,1);
if I==1, 
  if outs{1}(1)=='p',  
    fprintf(outs{1});
  end
else, 
  for i = 1:I,  
    if outs{i}{1}(1)=='p',
      J = size(outs{i},2);  
      for j = 1:J,  
        fprintf(outs{i}{j});
      end
    end
  end
end
end
%% print_report_fail
% [] = print_report_fail(outs)
% prints each line of output results
%

function [] = print_report_fail(outs,do_DNE,do_only_DNE), 
I = size(outs,1);
if I==1, 
  if outs{1}(1)=='F',  
    Ds = find(outs{1}(:)=='D');
    Ns = find(outs{1}(:)=='N');
    Es = find(outs{1}(:)=='E');
    if length(Ds) && length(Ns) && length(Es), 
      if do_DNE,  
        fprintf(outs{1});
      end 
    else, 
      fprintf(outs{1});
    end
  end
else, 
  for i = 1:I,  
    if outs{i}{1}(1)=='F',
      Ds = find(outs{i}{1}(:)=='D');
      Ns = find(outs{i}{1}(:)=='N');
      Es = find(outs{i}{1}(:)=='E');
      if length(Ds) && length(Ns) && length(Es), 
        if do_DNE,  
          J = size(outs{i},2);  
          for j = 1:J,     
              fprintf(outs{i}{j});
          end
        end
      else, 
        if ~do_only_DNE,  
          J = size(outs{i},2);  
          for j = 1:J,     
              fprintf(outs{i}{j});
          end
        end
      end
    end
  end
end
end
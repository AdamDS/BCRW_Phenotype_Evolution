%% generatlize_base_name.m
% [base_name] = generalize_base_name(base_name)
function [base_name] = generalize_base_name(base_name)
global SIMOPTS;
if SIMOPTS.exp_type==0, 
  us = 1;
  if SIMOPTS.reproduction==1, 
    us = 2;
  elseif SIMOPTS.reproduction==2, 
    us = 3;
  end
  if SIMOPTS.mu>=1.0, 
    ns = us +2;
  else
    ns = us +1;
  end
  loc_us = find(base_name=='_');
  bbn = base_name(1:loc_us(us));
  ebn = base_name(loc_us(ns):length(base_name)-1);
  base_name = [bbn 'x' ebn];
end
end
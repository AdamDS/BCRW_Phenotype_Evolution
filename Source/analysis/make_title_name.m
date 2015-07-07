%% make_title_name.m
% [title_name] = make_title_name(base_name,run_name)
%
% This function takes LaTeX into consideration for titles by replacing
% basic underscores with LaTeX underscores (i.e. _ to \_).
function [title_name] = make_title_name(base_name,run_name)
exp_name = [base_name run_name];
i = 1;  usi = 1;
title_name = exp_name;
us = find(title_name=='_');
while usi<=length(us)
  us = find(title_name=='_');
  okay = title_name(1:(us(usi)-1));  shift = title_name((us(usi)+1):length(title_name));
  title_name = [okay '\_' shift];
  usi = usi +1;
end
end
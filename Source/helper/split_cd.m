%% split_cd.m
% [new_dir_name] = split_cd(dir_name,run,split,make_dir,do_cd)
% Determines source directory based on split and run if data is hashed up
% according to simulation sorting in folders under the *\Data\ subtree.
% NOTE: Numbering starts at 0, although no run value may actually be 0.
% This should not cause problems.
%
% GENVARS = struct('babies',[],'basic_map',[],'land',[],'run',[],...
%   'exp_name',[],'dir_name',[]);
%
% GENOUTS = struct('population',population,'trace_noise',trace_noise,'trace_x',trace_x,...
%   'trace_y',trace_y,'trace_cluster_seed',trace_cluster_seed,'seed_distances',seed_distances,...
%   'parents',parents,'kills',kills,'rivalries',rivalries,'land',land,'shifted',shifted,...
%   'finished',finished);
%
function [new_dir_name] = split_cd(dir_name,run,split,make_dir,do_cd),
level = floor(run/split); % Hash up every 10 sims?
new_dir_name = dir_name;
if split, 
  new_dir_name = [dir_name ...
  int2str(level*split) '_' int2str(((1+level)*split)-1) '\'];
  if make_dir,  
    if exist(new_dir_name)==7,  
      if do_cd, 
        cd(new_dir_name);
      end
    else, 
      mkdir(new_dir_name);
      if do_cd, 
        cd(new_dir_name);
      end
    end
  else, 
    if exist(new_dir_name)==7,  
      if do_cd, 
        cd(new_dir_name);
      end
    end
  end
end
end
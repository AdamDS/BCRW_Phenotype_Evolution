%% get_gyration_radii.m
%
% Fixes gyration_radii, for a set of data. The orgsnclusters data is
% required to make the fix. Before 22 Aug 2013, gyration_radii was
% calculated as: 
% radii(cluster) = sqrt(mean((sum((these_coords-these_cents).^2,2))./orgs));
% when they should have been calculated as:
% radii(cluster) = sqrt(mean(sum((these_coords-these_cents).^2,2)));
%
% -ADS 10*4*12
function [gyration_radii] = fix_gyration_radii(base_name,run,dir_name), 
global SIMOPTS;
this_function = 'fix_gyration_radii';

gyration_radii = [];

new_dir_name = split_cd(dir_name,run,SIMOPTS.split,1,0);
run_name = int2str(run);
[clus_name] = cluster_name(base_name);
gr_name = [new_dir_name 'gyration_radii_' clus_name run_name];
onc_name = [new_dir_name 'orgsnclusters_' clus_name run_name];

fprintf(['Attempting ' this_function ' for ' clus_name run_name '\n']);


[gr,go,~] = try_catch_load([gr_name '.mat'],1,1);
if go,  [onc,go,~] = try_catch_load([onc_name '.mat'],1,1);
if go
  gyration_radii = gr.gyration_radii; clear gr
  orgsnclusters = onc.orgsnclusters';  clear onc

  gyration_radii = gyration_radii.*sqrt(orgsnclusters);
      %gyration radius for distances between center of cluster mass and
      %each organism (Lesne p15 assuming integral A(xbar) = orgs)
%       radii(cluster) = sqrt(mean((sum((these_coords-these_cents).^2,2))./orgs));
%FIXED ON 22 AUG 2013; ALL gyration_radii DATA BEFORE NEEDS TO BE FIXED!!!
%       radii(cluster) = sqrt(mean(sum((these_coords-these_cents).^2,2)));
  save([new_dir_name 'gyration_radii_' clus_name run_name],'gyration_radii');
%   save([new_dir_name 'gyration_radii2_' clus_name run_name],'gyration_radii2');
end %orgsnclusters
end %gyration_radii
end %function
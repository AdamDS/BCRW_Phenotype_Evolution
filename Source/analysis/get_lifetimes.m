%% get_lifetimes.m *******************************************************
% The output, ngens, is organized such that rows correspond to lifetimes results
% from a number of simulations, and the columns are mutablities grouped by
% death_max values. Differing sets of death_max values are separated by a
% column of zeros.

% I could set this up in an alternative way.
% Since the indiv_death_max simulations are run for many death_max values
% for each mutability, you can swap the typical loop order of mutability and
% random death values. This will produce the desired output, ngens, that is
% organized by simulation result in the rows and death_max value along the
% columns. Different mutability may be run simulataneously; however, the
% results from different mutabilities will be separated by a column of
% zeros. I would suggest that if you need many bn and death_max value sets, 
% then run each separately or in small groups.
global SIMOPTS;
ngens = []; %initialize the output
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
  i_ngens = [];
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  lt_name = make_data_name('lifetimes',base_name,'',0);
  load(lt_name); 
  lt = lifetimes; clear lifetimes
  lifetimes = lt; clear lt
  i_ngens = cat_row(i_ngens,lifetimes');
end
  if size(ngens,1)>0
    ngens = cat_row(ngens',zeros(size(i_ngens,2),1));
    ngens = cat_row(ngens,i_ngens'); %
  else
    ngens = i_ngens';
  end
end
end
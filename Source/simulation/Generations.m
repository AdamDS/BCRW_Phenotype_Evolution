%% Generations.m
% function [GENOUTS] = Generations(GENVARS)
% inputs:
% GENVARS {babies, land, basic_map, exp_name}
%
% outputs:
% GENOUTS {population, trace_noise, trace_x, trace_y, trace_cluster_seed, seed_distances,...
%          parents, kills, rivalries, land, shifted, finished}
%
%   mnd = zeros(SIMOPTS.NGEN,1); %records the mean neighbor distance in each generation
%   mnd = [];
% trace_cluster = []; %records cluster number of traced animals
%
function [GENOUTS] = Generations(GENVARS)
global SIMOPTS;
%% Initialize 
  % Output data
population = zeros(1,SIMOPTS.NGEN); %records the population of each generation
if SIMOPTS.save_kills, 
  kills = zeros(SIMOPTS.NGEN,3);
  if SIMOPTS.save_rivalries,  
    rivalries = zeros(SIMOPTS.NGEN,1); %records the number of sibling rivalries occuring in each generation
  else, rivalries = []; end
else, %if only_lt, then there is no need to record these; this leaves them empty
  kills = []; %records the kill counts from each death type
  rivalries = [];
end
trace_cluster_seed = []; %records mate and alternate of each organism
seed_distances = []; %records the nearest & second nearest neighbor distances to each organism
trace_x = []; %records x position of traced animals
trace_y = []; %records y position of traced animals
trace_noise = []; % records noise value of traced animals
shifted = []; %records shifted landscape
parents = []; %records parent(s) of each organism
adults = []; %generational record of parents
finished = SIMOPTS.NGEN; %if simulation ends early, this will tell the final generation simulated
babies = GENVARS.babies;
land = GENVARS.land;
basic_map = GENVARS.basic_map;

% loop parameters
gen = 0; %track generation
done = 0; %determines when loop is done
while ~done, 
  gen = gen +1; %increment the generation
    
%% Report simulation progress to the Command Window
  print_simulation_update(gen,size(babies,1,...
    numel(unique(babies(:,3))),SIMOPTS.exp_name,SIMOPTS.run);

%% Determine nearest (mates) and second nearest (alternates) neighbors
  [M,M_DIST,SN,SN_DIST] = FindMates(babies,land);

%% Record basic population data
  population(gen) = size(babies,1); %record population
  if ~SIMOPTS.only_lt, 
    if SIMOPTS.record, 
      if SIMOPTS.exp_type~=0, trace_noise = [trace_noise; babies(:,3)]; end;  %record mutability
      trace_x = [trace_x; babies(:,1)]; %records x position of traced animals
      trace_y = [trace_y; babies(:,2)]; %records y position of traced animals
    end
    if SIMOPTS.save_seeds,  
      trace_cluster_seed = [trace_cluster_seed; [M,SN]]; %record cluster seeds
      seed_distances = [seed_distances; [M_DIST,SN_DIST]];
    end
    if SIMOPTS.save_parents,  parents = [parents; adults];  end
  end

%% Reproduction & Death
  %BABIES GROW UP
  par = babies; %current babies are the new parents!

  %MAKE BABIES
  babies=[]; %clear previous population of babies
  [babies,adults] = MakeBabies(par,M,land);

  %OVERPOPULATION DEATH
  [babies,adults,odkills,sibrivals] = OverpopulationDeath(babies,adults,land);

  %ALSO SOME RANDOM DEATH FOR THIS GENERATION
  [babies,adults,rdkills] = RandomDeath(babies,adults);

  %ALSO KILL BABIES OUTSIDE OF phenotype values allowed in the mapping
  if SIMOPTS.reproduction~=3,  
    [babies,adults,cjkills] = CliffJumpers(babies,adults,land);
  else, cjkills = []; end

  %ADJUST LANDSCAPE randomly for the babies.
  [land,basic_map,old_land] = AdjustLandscape(basic_map,land,gen);
  %%% Does not have the periodic boundaries update. %%%
  %%% Does not have feedback. %%%
    
  if ~SIMOPTS.only_lt, 
    if SIMOPTS.save_kills,  
      %Collect kill counts
      kills(gen,:) = [odkills, rdkills, cjkills];
      if SIMOPTS.save_rivalries,  rivalries(gen) = [sibrivals]; end
    end
    %Collect landscapes
    if ~SIMOPTS.shock,  shifted = [old_land; shifted(2:(size(shifted,1)),:)]; %update shifted
    else, shifted = [old_land; shifted]; %update shifted
    end
  end
%% Determine whether to end generation looping
  %if competing mus & only 1 mu left      or last gen  or not enough babies     then end of sim                   
  if ((SIMOPTS.exp_type==1 && length(unique(babies(:,3)))==1) ...
     || (gen==SIMOPTS.NGEN) || (size(babies,1)<SIMOPTS.limit)), 
    if (SIMOPTS.IPOP<SIMOPTS.limit && size(babies,1)<=0) || ~SIMOPTS.loaded,  
      done = 1; 
    end
  end
end
finished = gen; %report the last generation simulated
%% Clean up some output data
if ~SIMOPTS.only_lt, 
  population = population(population~=0 | population(gen)>=SIMOPTS.limit); %get non-zero populations
  rivalries = rivalries(population>=SIMOPTS.limit);
  %In case simulation ends earlier than NGEN
  if finished<SIMOPTS.NGEN, 
    %get only non-zero rows of shifted
    if SIMOPTS.landscape_heights(1)~=SIMOPTS.landscape_heights(2),  
      shifted = shifted(shifted(:,1)~=0);
    end
  end
end

%% Bundle the output data into GENOUTS
GENOUTS = struct('population',population,'trace_noise',trace_noise,'trace_x',trace_x,...
  'trace_y',trace_y,'trace_cluster_seed',trace_cluster_seed,'seed_distances',seed_distances,...
  'parents',parents,'kills',kills,'rivalries',rivalries,'land',land,'shifted',shifted,...
  'finished',finished);
end
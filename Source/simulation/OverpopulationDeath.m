% OverpopulationDeath.m
% function [babies,parents] = OverpopulationDeath(babies,overpop,parents)
% inputs:
% babies - offspring to check
% parents - parents of the offspring being checked
% land - landscape, used with periodic boundaries only
%
% outputs:
% babies - offspring which survived
% parents - parents of the offspring that survived
% odkills - overpopulation death kill count
% sibrivals - instances of sibling rivalry
%
% This function has been updated from its orginal version to:
%  -randomize the order of reference offspring to be checked for overpopulation death
%  -include save options for parents and kill count
%  -include an option to choose from random walk styles: coalescing (original) or annihilating
%  -determine distances in the case of periodic boundaries
% Original: function [babies] = OverpopulationDeath(babies,overpop)
% function [babies,parents,odkills,sibrivals] = OverpopulationDeath(babies,overpop,parents,...
%           save_parents,save_kills,random_walk,save_rivalries)
function [babies,parents,odkills,sibrivals] = OverpopulationDeath(babies,parents,land)
global SIMOPTS;
%% Initialize a few parameters
op2 = SIMOPTS.op^2; %since no sqrt's to save computation time
mark = 10^4; %a marker which will certainly remove organisms from the landscape

%% Randomize the order of offspring which act as reference organisms
rand_baby_list = rand(size(babies,1),1); %assign each baby a random number
[~, I] = sort(rand_baby_list); %srlist is the ordered list of random assignments
                       %I is the ordered list of indices of rand_baby_list
kids = babies(I,:); %randomized order of offspring
sibrivals = 0;

wrap_up = size(land) +0.5 -SIMOPTS.op;  %furthest wrappables
wrap_down = 0.5+SIMOPTS.op; %nearest wrappables
  %the +0.5 is because the organisms can exist between 0.5 & size(land)+0.5
  %extending inward by an amount of SIMOPTS.op accounts for the boundary of
  %potential crowders from across the boundaries
mod_d2 = size(land).^2 -op2;
if SIMOPTS.periodic(1), wrappers_x = find(kids(:,1)>wrap_up(1) | kids(:,1)<wrap_down);
else, wrappers_x = [];  end
if SIMOPTS.periodic(2), wrappers_y = find(kids(:,2)>wrap_up(2) | kids(:,2)<wrap_down);
else, wrappers_y = [];  end

if SIMOPTS.random_walk==0,  
  %% Coalescing
  for i = 1:size(babies,1), %for each baby
    if kids(i,1)~=mark, %start with a live one
      not_dead = find(kids(:,1)~=mark); %don't work on dead guys
      not_self = find(i~=not_dead); %don't calculate distance to self
      dx2 = (kids(not_dead(not_self),1) -kids(i,1)).^2;
      dy2 = (kids(not_dead(not_self),2) -kids(i,2)).^2;
      if find(i==wrappers_x), dx2 = mod(dx2,mod_d2(1));  end
      if find(i==wrappers_y), dy2 = mod(dy2,mod_d2(2));  end
      baby_distance = dx2 + dy2; %find distance to all babies
      %%% could determine a small range instead of checking the distances to every other baby %%%
      [kill]=find(baby_distance < op2); %find those that are too close
%to track sibling rivalry
      [~, dis_guy] = min(baby_distance(baby_distance>=op2)); %get nearest okay neighbor
%       dnn(not_dead(not_self(dis_guy(1)))) = distnn(1); %trace nearest okay distance
      if SIMOPTS.save_rivalries && ~isempty(kill), %if recording sibling rivalries & killing someone
        par_of_ref_off = parents(I(i),:); %parents of reference
        par_of_kill_off = parents(I(not_dead(not_self(kill))),:); %parents of killings
        rival_siblings = find(par_of_ref_off(1,1)==par_of_kill_off(:,1)); %determine if siblings
        if ~isempty(rival_siblings), %if they are
          sibrivals = sibrivals +1; %increment sibling rivalry kills
%         else
%           kill = []; %if not sibling rivalry, then no death (use for control conditions only)
        end
        rival_siblings = []; %refresh
      end
%end tracking of sibling rivalries
      %mark those too close, not_dead(not_self(kill))
      kids(not_dead(not_self(kill)),:) = mark*ones(length(kill),3);
    end
  end %end of overpopulation death
elseif SIMOPTS.random_walk==1,  
  %% Annihilating
  for i = 1:size(babies,1), %for each baby that is still alive (number is changing...)
    if kids(i,1)~=mark, %don't worry about those already marked for death
      not_dead = find(kids(:,1)~=mark); %find those not yet marked for death
      not_self = find(i~=not_dead); %don't worry about self distance (which is always 0)
      dx2 = (kids(not_dead(not_self),1) -kids(i,1)).^2;
      dy2 = (kids(not_dead(not_self),2) -kids(i,2)).^2;
      if find(i==wrappers_x), dx2 = mod(dx2,mod_d2(1));  end
      if find(i==wrappers_y), dy2 = mod(dy2,mod_d2(2));  end
      baby_distance = dx2 + dy2; %find distance to all babies
      %%% could determine a small range instead of checking the distances to every other baby %%%
      [kill]=find(baby_distance < op2); %find those that are too close
%to track sibling rivalry
      [~, dis_guy] = min(baby_distance(baby_distance>=op2)); %get nearest okay neighbor
%       dnn(not_dead(not_self(dis_guy(1)))) = distnn(1); %trace nearest okay distance
      if SIMOPTS.save_rivalries && ~isempty(kill), %if recording sibling rivalries & killing someone
        par_of_ref_off = parents(I(i),:); %parents of reference
        par_of_kill_off = parents(I(not_dead(not_self(kill))),:); %parents of killings
        rival_siblings = find(par_of_ref_off(1,1)==par_of_kill_off(:,1)); %determine if siblings
        if ~isempty(rival_siblings), %if they are
          sibrivals = sibrivals +1; %increment sibling rivalry kills
%         else
%           kill = []; %if not sibling rivalry, then no death (use for control conditions only)
        end
        rival_siblings = []; %refresh
      end
%end tracking of sibling rivalries
      if ~isempty(kill),  
        %mark the reference organism, i, and those too close, not_dead(not_self(kill))
        kids([i; not_dead(not_self(kill))],:) = mark*ones((length(kill) +1),3);
      end
    end
  end %end of overpopulation death
end
if SIMOPTS.save_kills,  
  odkills = length(find(kids(:,1)==mark)); 
else, 
  odkills = [];  
end %record kill count
babies(I,:) = kids; %update babies with the de-randomized list (sort(I))
% nn_dist = dnn(find(dnn(sort(I))));
kill_them = find(babies(:,1)==mark); %determine who has been marked for death
babies(kill_them,:) = []; %kill those marked for death
if SIMOPTS.save_parents,  parents(kill_them,:) = [];  end %remove the dead organisms' parents
end
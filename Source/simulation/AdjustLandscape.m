%% AdjustLandscape.m
% inputs - basic_map, landscape
% This has been updated to do Shocks
% This has been updated to work with global SIMOPTS.
% This is NOT updated for periodic boundaries.
% function [land,basic_map,old_land] = AdjustLandscape(basic_map,land,gen)
function [land,basic_map,old_land] = AdjustLandscape(basic_map,land,gen)
global SIMOPTS;
old_land = [];
if ~SIMOPTS.shock, %if no shock
  if ~mod(gen,SIMOPTS.landscape_movement),
    if SIMOPTS.landscape_heights(1)~=SIMOPTS.landscape_heights(2), %uneven landscape
      old_land = land(end-4:end,:);
      height_difference = SIMOPTS.landscape_heights(2) -SIMOPTS.landscape_heights(1);
      basic_map=[rand(1,size(basic_map,1))*height_difference +SIMOPTS.landscape_heights(1); 
      basic_map(1:end-1,:);]; %create new random landscape column for left side
      land=ceil(interp2(basic_map,2)); %put new map together and interpolate
    end
  end %end of landscape adjustment
else, %if shocking
  if ~SIMOPTS.shock_over, 
    if SIMOPTS.shock_duration~=Inf, %finite shock
      if ~SIMOPTS.shock_repeat, %shock once
        if gen>=SIMOPTS.landscape_movement, %after beginning of shock
          if gen<SIMOPTS.landscape_movement +SIMOPTS.shock_duration, %during shock
            if SIMOPTS.shock_heights(1)==SIMOPTS.shock_heights(2), 
              basic_map = SIMOPTS.shock_heights(1)*ones(size(basic_map));
              land = SIMOPTS.shock_heights(1)*ones(size(land));
            else, 
              %no sims for uneven shock landscape yet
            end
          else, %after end of shock
            if SIMOPTS.landscape_heights(1)==SIMOPTS.landscape_heights(2), 
              SIMOPTS.shock_over = 1; %set flag to stop shocking
              basic_map = SIMOPTS.landscape_heights(1)*ones(size(basic_map));
              land = SIMOPTS.landscape_heights(1)*ones(size(land));
            else, 
              %no sims for uneven regular landscape yet
            end
          end
        end
      else, %repeated shocks
        %no sims for this yet
      end
    else, %forever shocking
      if ~mod(gen,SIMOPTS.landscape_movement),
        old_land = basic_map;
        basic_map=rand(size(basic_map))*SIMOPTS.landscape_heights(2);
        land=ceil(interp2(basic_map,2)); %put new map together and interpolate
      end
    end
  end
end
end

%
%   llllllllllll    nnnnnnnnnnnn    new row of basic_map
%   llllllllllll    ssssssssssss    of previous basic_map
%   llllllllllll    ssssssssssss
%   llllllllllll    ssssssssssss
%   llllllllllll    ssssssssssss
%   llllllllllll    ssssssssssss
%   llllllllllll    ssssssssssss
%   llllllllllll    ssssssssssss
%   llllllllllll    ssssssssssss
%   llllllllllll    ssssssssssss
%   llllllllllll    ssssssssssss
%   llllllllllll    ssssssssssss
%   llllllllllll    oooooooooooo    old row of basic_map
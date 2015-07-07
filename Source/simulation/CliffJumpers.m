%% CliffJumpers
% input - babies, parents, landscape
% output - babies that didn't go out of bounds, their parents, number killed
% Those babies which are not between 0.5 and max(size(land))+0.5 along either phenotype
% axis will be killed.  This code is taken directly from Nate's design
% along with the save parents option. 
% This has been updated to work with global SIMOPTS.
% This is updated for periodic boundaries.
% function [babies,adults,cjkills] = CliffJumpers(babies,adults,land)
function [babies,adults,cjkills] = CliffJumpers(babies,adults,land),  
global SIMOPTS;
outside_x = []; outside_y = []; outside = [];
if SIMOPTS.periodic(1), %if x is wrapped
  if SIMOPTS.periodic(2), %if y is wrapped
    outside = []; %no outsiders
  else, %find outside y
    outside = find(babies(:,2) > size(land,2)+0.5 | babies(:,2) < 0.5);
  end
else, %x isn't wrapped
  if SIMOPTS.periodic(2), %if y is wrapped
    outside = find(babies(:,1) > size(land,1)+0.5 | babies(:,1) < 0.5);
  else, %neither wrapped
    outside=find(babies(:,1) > size(land,1)+0.5 | babies(:,1) < 0.5 |...
                 babies(:,2) > size(land,2)+0.5 | babies(:,2) < 0.5);
  end
end
babies(outside,:) = []; %kill them.
if SIMOPTS.save_parents,  adults(outside,:) = []; end %remove their parents
if SIMOPTS.save_kills,   cjkills = length(outside); else, cjkills = []; end
end
%% Cluster edge finding
% This algorithm determines a boundary about a known set of points within a
% cluster. The only input needed are the x & y coordinates of the cluster
% points. The algorithm starts by finding the lowest y point. If multiple
% points are along the "bottom" then the largest x point of those is taken.
% This first point is the "seed". From here, the slopes connecting each
% other point to the seed is determined. Taking the minimum slope will give
% the next point. If multiple points are in line, then the most distant
% point is taken. The cluster points are all then translated an amount such
% that the next point sits at the origin. The cluster is then rotated by an
% amount such that the previous point sits along the negative x-axis. A
% search is done for any points to the right of the shifted cluster origin.
% This step is inserted to help speed up the algorithm. If no points are 
% found, then all points left are then considered. The points considered
% have their slopes to the origin measured, and the point making the lowest
% slope is taken. The algorithm repeats until the next point is the seed.
% The output is the list of indices of COORDS, "points" which lists the
% vertices of the polygon boundary & the list of the coordinates of
% "points", "edge".
% -ADS 2 December 2012
clear all; %close all;
do_plot = 1;  do_steps = 0;
N = 5;  %number of points 
lims = 1.5; %limits of the graph xlim([-lims lims]); ylim([-lims lims]);
% Comment out whichever set you don't want to work with.
COORDS = rand(N,2); %work with just random points?
% COORDS = [0 0; 1 0; 1 1; 0 1; 0.5 0.9]; %work with some known shape?
% load('broken');
coords = COORDS;
%% Begin algorithm!
% for asdf=1:1000, 
% tic;
N = size(coords,1);
edge = [];  e = [];
points = [];   p = [];
%% start with a certain point, minimum y
from = find(min(coords(:,2))==coords(:,2));
if numel(from)>1, % in case of in-line points
  [X from] = max(coords(from,1)); % get the far right point
  Y = coords(from,2);
  seed = from; %initialize seed
  %when completing the polygon, the far right point will be the one chosen
  %according to the algorithm, and since the loop checks if the points,
  %next & seed are the same, the loop should certainly exit with this
  %assignment.
else, 
  X = coords(from,1);
  Y = coords(from,2);
  seed = from; %initialize seed
end
next = -1;
%% main algorithm loop
while next~=seed, 
  % no points are below from, so check to the right for counter-clockwise
  % hull path
  to = find(X<=coords(:,1));
  to(to==from) = []; %don't include from
  if isempty(to), %if no points left, start with a full non-visited subset
    to = 1:N; %all points are left of current from coord
    to(from) = []; %don't include from
    theta = pi/2;
    big_bends = find(Y<=coords(to,2));
    if numel(to)==numel(big_bends), 
      theta = pi;
    end
  else, 
    theta = 0;
  end
  dy = (coords(to,2)-Y);
  dx = (coords(to,1)-X);
  slopes = (dy./dx);
  j = (find(min(slopes)==slopes));
  next = to(j);
  if numel(next)>1, % in case of in-line points
    next = max(coords(next,1));
  end
  smallest_slope = slopes(j);
if do_plot==1, 
  if ~figure(1),  figure(1);  end
  subplot(3,1,1);
  plot(COORDS(:,1),COORDS(:,2),'x'); hold on; plot(COORDS([seed p' from next],1),COORDS([seed p' from next],2)); hold off;
  xlim([-lims lims]); ylim([-lims lims]);
end
  translated_coords(:,1) = coords(:,1)-coords(next,1);
  translated_coords(:,2) = coords(:,2)-coords(next,2);
if do_plot==1, 
  subplot(3,1,2);
  plot(translated_coords(:,1),translated_coords(:,2),'x'); 
  hold on; plot(translated_coords([seed p' from next],1),translated_coords([seed p' from next],2)); hold off;
  xlim([-3 3]); ylim([-3 3]);
end
  if abs(dx(j))>10^-15, 
    theta = theta +atan(smallest_slope);
  else, 
    theta = theta +pi/2;
  end
  rotation_matrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
  coords = translated_coords*rotation_matrix;
if do_plot==1, 
  subplot(3,1,3);
  plot(coords(:,1),coords(:,2),'x');
  hold on; plot(coords([seed p' from next],1),coords([seed p' from next],2)); hold off;
  xlim([-lims lims]); ylim([-lims lims]);
end
if do_steps==1, pause;  end
  X = coords(next,1);
  Y = coords(next,2);
  % record point id and coordinate
  if from~=seed, 
    p = [p; from];
    e = [e; COORDS(from,:)];
  end
  % update next point
  from = next;  
end
points = [seed; p; seed];
edge = [COORDS(seed,:); e; COORDS(seed,:)];
% times(asdf) = toc;
% end
[points, edge]
if ~figure(2),  figure(2);  end
plot(edge(:,1),edge(:,2)); xlim([-lims lims]); 
ylim([-lims lims]);
hold on;
plot(COORDS(:,1),COORDS(:,2),'x');
hold off;
%% old notes to myself -ADS
% I want this algorithm to determine an edge of a predetermined cluster.
% Given a simple shape, the algorithm should correctly determine the shape.
% I will first consider a counter-clockwise path about a cluster, but
% perhaps later it can be modified to trace clockwise as well. I think I
% will first try a scan of points from the outer section inward. Once a
% point is found which is a potential next vertex of the polygon being
% traced, I want any points in between and within a circle of radius equal
% to the distance from the current vertex to the next potential vertex. Of
% these points in between, which should only contain points "more inward"
% (if in the end they/it is included, then they/it would form what might be
% considered a dent), if the distance to them is less than the distance to
% the potential vertex, then they become a part of the edge. This should
% cause a tight fitting polygon which should be able to pick up on dents.
%
% Instead of the above, I think a more tractable algorithm could be
% obtained by use of Euler rotations and a greatest slope search. First, a
% search is done for the next vertex which makes the lowest slope of a line
% stemming form the current vertex, but only points beyond the current vertex
% are searched. The cluster of points is then rotated by an angle such that 
% the slope is reduced to zero.
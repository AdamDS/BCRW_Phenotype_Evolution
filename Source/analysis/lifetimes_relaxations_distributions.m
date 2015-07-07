clear all; clc; close all;
dmu = 0.001; mutabilities = [0.332:dmu:0.332];
mult = 1/dmu;
NGEN = 10^5;

IPOP = 30000; % 6x6=6549, 8x8=12474, 10x10=20290, 12x12=30000, 20x20=87757

basic_map_size = [12 12];

ngen = [int2str(NGEN) '_NGEN_'];
ipop = [int2str(IPOP) '_IPOP_'];
if basic_map_size(1)==12 && basic_map_size(2)==12,  bms = ''; 
else, bms = [int2str(basic_map_size(1)) 'x' int2str(basic_map_size(2)) '_basic_map_size_']; end

ds = 18;  fin = 100;
start = [1:ds:fin];
last = [ds:ds:fin 100];

for m = mutabilities, 
  fprintf('Combining mu = %1.4f\n',m);
  if mod(m*mult,10),  
    mu = int2str(m*mult);  %subplot(3,2,4);
  else, 
    mu = int2str(m*10);
  end
  % s0 = 330; split = 5;  
  do_lt_save = 0;
  do_p_save = 0;
  lt = [];  p = [];
  avail = 1:length(start);
  for i = avail, 
  %   start = s0+split*(i-1);  last = (s0-1)+split*i;
    s = int2str(start(i)); l = int2str(last(i)); 
    range = [s '_' l];
    if do_lt_save==(i-1) || i==1,  
      [lifes,have,~] = try_catch_load(['lifetimes_Bacterial_Mu_' mu '_Uniform_2_Flatscape_' ...
        bms ipop ngen range '.mat'],1,1);
      if have,  lt = [lt; lifes.lifetimes]; do_lt_save = 1 +do_lt_save; end
    end
    if do_p_save==(i-1) || i==1,  
      [pops,have,~] = try_catch_load(['populations_Bacterial_Mu_' mu '_Uniform_2_Flatscape_' ...
        bms ipop ngen range '.mat'],1,1);
      if have,  p = cat_row(p,pops.populations); do_p_save = 1 +do_p_save;  end
    end
  end

  if do_lt_save==length(avail), 
    beg = start(1);  las = l;
    ran = [int2str(beg) '_' l];
    lifetimes = lt;
    save(['lifetimes_Bacterial_Mu_' mu '_Uniform_2_Flatscape_' bms ipop ngen ran '.mat'],'lifetimes');
  end
  if do_p_save==length(avail),  
    beg = start(1);  las = l;
    ran = [int2str(beg) '_' l];
    populations = p;
    save(['populations_Bacterial_Mu_' mu '_Uniform_2_Flatscape_' bms ipop ngen ran '.mat'],'populations','-v7.3');
  end

end
% figure, hist(lifetimes);
% figure, plot(mean(populations,1));
% title(['\mu = 0.' mu],'FontSize',20);
% xlabel('generations','FontSize',20);
% ylabel('frequency','FontSize',20);

% mu = int2str(325);  subplot(2,2,4);
% load(['relaxations_Bacterial_Mu_' mu '_Uniform_2_Flatscape_30000_IPOP_10000000_NGEN.mat']);
% hist(relaxations);
% title(['\mu = 0.' mu],'FontSize',20);
% xlabel('generations','FontSize',20);
% ylabel('frequency','FontSize',20);
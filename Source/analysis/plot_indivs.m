%plot_indivs.m
% This function allows you to specify a string identifying a simulation set, a string
% telling which specific simulation, and for what generation you'd like to see where
% organisms are in the phenospace.
% Inputs:
% base_name = string identifying a simulation set
% run_name = string identifying the particular simulation number
% gen = generation of organisms to plot
% function [] = plot_indivs(base_name,run_name,gen)
% function [] = plot_indivs(base_name,run_name,gen,color_code)
global SIMOPTS;
m = 0;  cc = color_code;  ngens = [];
if backgroundcolor=='k',  markercolor = 'w';  
elseif backgroundcolor=='w',  markercolor = 'k';  end
if generate_paper_plot || prod(subplotrc)==9, h = 0;  end
shapes = ['.+sv'];  %+o*.xsdv^><ph

for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  m = m +1;
  if mu==plot_mu(m),  
  for run = SIMS, 
    run_name = int2str(run);
    color_code = cc;
    if run==plot_run
    if generate_paper_plot || prod(subplotrc)==9, 
      fign = ceil(1000*mutability(1) +10^6*mutability(2) +10^9*mutability(3));
    else, 
      if limit~=3
        fign = ceil(1000*mu)+limit;
      else
        fign = ceil(1000*mu);
      end
    end
    figure(fign);
    if prod(subplotrc)==4,  set(figure(fign),'WindowStyle','normal','Position',[1 50 950 950]);
    elseif prod(subplotrc)==3,  
      set(figure(fign),'WindowStyle','normal','Position',[200 200 810 270]);  
    elseif prod(subplotrc)==9,  
      set(figure(fign),'WindowStyle','normal','Position',[200 200 365 365]);  
    end
    set(figure(fign),'WindowStyle','normal','Position',[00 00 1000 1000]); 
%     if ~generate_paper_plot || prod(subplotrc)~=9, h = 0;  end
    go = 1;
    [p,go,~] = try_catch_load(['population_' base_name run_name],go,1);
    if go,  population = p.population;  clear p
    [tx,go,~] = try_catch_load(['trace_x_' base_name run_name],go,1);
    if go,  trace_x = tx.trace_x; clear tx
    [ty,go,~] = try_catch_load(['trace_y_' base_name run_name],go,1);
    if go,  trace_y = ty.trace_y; clear ty
    [tc,gotc,~] = try_catch_load(['trace_cluster_' base_name run_name],go,1);
    if gotc,  trace_cluster = tc.trace_cluster; clear tc, 
      [onc,goonc,~] = try_catch_load(['orgsnclusters_' base_name run_name],go,1);
      if goonc, orgsnclusters = onc.orgsnclusters;  clear onc,  end
      [nc,gonc,~] = try_catch_load(['num_clusters_' base_name run_name],go,1);
      if gonc,  num_clusters = nc.num_clusters; clear nc, end
    end
    ngen = length(population);
    ngens = [ngens ngen];
    if max(plot_gens)>ngen, plot_gens = [1 ceil(ngen/2) ngen];  end
    if ~subplots && make_video
      plot_gens = 1:ngen;
    end
    for generation = plot_gens, 
      %generate new subplot space
      if subplots,  h = h +1; subplot(subplotrc(1),subplotrc(2),h); end

      v = sum(population(1:generation));  u = v -population(generation) +1;
%       t = [];
%       if gotc && gonc && goonc,  
%         if plot_largest_clusters, 
%           t = u:v;
%           cv = sum(num_clusters(1:generation)); cu = cv -num_clusters(generation) +1;
%           [~,isonc] = sort(orgsnclusters(cu:cv));
%           T = [];
%           nctp = length(color_code);  %Number of Clusters To Plot
%           if nctp<max(trace_cluster(u:v), nctp = max(trace_cluster(u:v)); end
%           for j = length(color_code), 
%             t = [t find(trace_cluster(u:v)==isonc(end -j))'];
%           end
%         else, 
%           for cci = 1:length(color_code), 
%             t = [t find(trace_cluster(u:v)~=color_code(cci))'];
%           end
%         end
%         t = unique(t)' +u -1;
%       else, t = u:v;  end

%       plot(trace_x(t)-0.5,trace_y(t)-0.5,['.' markercolor],'MarkerSize',4);  hold on;
%       %plot(trace_x(t)-0.5,trace_y(t)-0.5,'.k','MarkerSize',4);  hold on;
%       xywh = get(gca,'Position');
%       set(gca,'Color',backgroundcolor,'PlotBoxAspectRatioMode','manual');
      if prod(subplotrc)==4,  
        if h==1,  set(gca,'Position',[0.010 0.510 0.45 0.45]);
        elseif h==2,  set(gca,'Position',[0.510 0.510 0.45 0.45]);
        elseif h==3,  set(gca,'Position',[0.010 0.010 0.45 0.45]);
        elseif h==4,  set(gca,'Position',[0.510 0.010 0.45 0.45]);  end
      elseif prod(subplotrc)==3,  
        if h==1,  set(gca,'Position',[0.0 0.0 0.3333 1.0]);
        elseif h==2,  set(gca,'Position',[0.333325 0.0 0.3333 1.0]);
        elseif h==3,  set(gca,'Position',[0.66665 0.0 0.3333 1.0]); end
      elseif prod(subplotrc)==9,  
        if h==7,  set(gca,'Position',[0.0 0.0 0.3333 0.3333]);
        elseif h==8,  set(gca,'Position',[0.333325 0.0 0.3333 0.3333]);
        elseif h==9,  set(gca,'Position',[0.66665 0.0 0.3333 0.3333]);
        elseif h==4,  set(gca,'Position',[0.0 0.333325 0.3333 0.3333]);
        elseif h==5,  set(gca,'Position',[0.333325 0.333325 0.3333 0.3333]);
        elseif h==6,  set(gca,'Position',[0.66665 0.333325 0.3333 0.3333]);
        elseif h==1,  set(gca,'Position',[0.0 0.66665 0.3333 0.3333]);
        elseif h==2,  set(gca,'Position',[0.333325 0.66665 0.3333 0.3333]);
        elseif h==3,  set(gca,'Position',[0.66665 0.66665 0.3333 0.3333]); end
      end

      %% Plot general population
%         colors = ['oworobdwdrdb'];
      plot(trace_x(u:v)-0.5,trace_y(u:v)-0.5,[shapes(1) markercolor],'MarkerSize',4);  hold on;
      set(gca,'Color',backgroundcolor);

      %% Plot specific clusters
      if generate_paper_plot, colors = ['skdkxkok'];  
      else, colors = ['oroboksrsbsk'];  faces = 'rbyrby'; end
      if ~generate_paper_plot || prod(subplotrc)~=9, 
        if max(trace_cluster(u:v))<max(color_code), %have enough clusters for available color_codes
          if length(color_code)>1, %if there are color_codes
            if max(trace_cluster(u:v))>1, %if there is more than one cluster
              l = 0;
              for k = 1:(color_code(2)-color_code(1)):color_code(max(trace_cluster(u:v))), 
                l = l +1;
                color_code(l) = k;
              end
            else, 
              cc = color_code;
              color_code = 1;
            end
          end
        end
        if length(color_code)>6, 
          color_code(7:length(color_code)) = [];
        end
        if color_code(1)>0, 
          j = 0;  fc = 0;
          for i = color_code, 
            j = j +2; fc = fc +1;
            ci = find(trace_cluster(u:v)==i);
            plot(trace_x(u-1+ci)-0.5,trace_y(u-1+ci)-0.5,colors(j-1:j),'MarkerSize',4,...
              'MarkerFaceColor',faces(fc));%,...
%                  'MarkerFaceColor',colors(j));
          end
        end
        
      else, 
        if plot_largest_clusters,
          cv = sum(num_clusters(1:generation)); cu = cv -num_clusters(generation) +1;
          [~,isonc] = sort(orgsnclusters(cu:cv));
          nctp = length(color_code);  %Number of Clusters To Plot
          if nctp>max(trace_cluster(u:v)), nctp = max(trace_cluster(u:v)); end
          for i = 0:(nctp-1), 
            these = find(trace_cluster(u:v)==isonc(end -i)) +u -1;
%             if shapes(i+2)=='s',  msize =3;  else, msize = 3;  end
            msize = 2;
            plot(trace_x(these)-0.5,trace_y(these)-0.5,[shapes(i+2) markercolor],'MarkerSize',msize,...
                 'MarkerFaceColor',markercolor);
          end
        end       
      end

      %% Format plots
      hold off;
      land_size = ((((basic_map_size*2)-1)*2)-1);
      buffer = 0.5;
      xlim([-buffer land_size(1)+buffer]); ylim([-buffer land_size(2)+buffer]);
      [tn] = make_title_name(base_name,run_name);

      if ~generate_paper_plot,  
%         title([make_title_name(base_name,run_name) ' @ generation ' int2str(generation)]);  
        % tick mark only on bottom left subplot figure
        if prod(subplotrc)==4, 
          if m==3 && i, 
            set(gca,'XTick',[0:floor(land_size(1)/3):land_size(1)],'XTickLabel',{},...
              'YTick',[0:floor(land_size(2)/3):land_size(2)],'YTickLabel',{});
          else, 
%             set(gca,'XTick',[0 land_size(1)],'XTickLabel',{},...
%               'YTick',[0 land_size(2)],'YTickLabel',{});
          end
        end
      else, %clean up axes and box plots
        set(gca,'XTick',[],'YTick',[],'Box','on');
      end
  set(gca,'XTick',[],'YTick',[],'Box','on');
      %% Make videos if desired
      if make_video, indiv_movie(generation) = getframe; end
    end
    end
    end
    end
    end
  end
    if make_video, movie2avi(indiv_movie,[base_name int2str(plot_run)],'compression','none'); end
  end
end
end
end
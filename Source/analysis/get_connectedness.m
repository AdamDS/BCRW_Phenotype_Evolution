global SIMOPTS;
sml = [];
sal = [];
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  [base_name,dir_name] = NameAndCD(make_dir,do_cd);
  clus_name = cluster_name(base_name);
  for run = SIMS
    %% Load data
    run_name = int2str(run);
    new_dir_name = split_cd(dir_name,run,split,make_dir,do_cd);
    [p,go] = try_catch_load([new_dir_name 'population_' base_name run_name],1,1);
    if go,  population = p.population;  clear p,  
    [tcs,go] = try_catch_load([new_dir_name 'trace_cluster_seed_' clus_name run_name],1,1);
    if go,  trace_cluster_seed = tcs.trace_cluster_seed;  clear tcs,  
      %% Initialization
      ngen = length(find(population>=limit));
  %     mate_links = zeros(ngen,10000,'uint16');
  %     if size(trace_cluster_seed,2)==2, second_links = zeros(ngen,10000,'uint16'); 
      mate_links = zeros(ngen,max(population));
      if size(trace_cluster_seed,2)>1, second_links = zeros(ngen,max(population)); 
      else second_links = []; end
      %% main algorithm
      for gen = 1:ngen, 
        [u,v] = gen_index(population,gen);
%         if population(gen)>=limit,  %make sure population is sufficient
        mates = trace_cluster_seed(u:v,1);  %get mates
        hm = hist(mates,1:length(mates)); %mate degree connectivity
        nzhm = find(hm~=0); %find non-zero frequencies
        mate_links(gen,nzhm) = hist(hm(nzhm),1:length(nzhm)); %frequency of mate degree
        if size(trace_cluster_seed,2)>1,  
          second = trace_cluster_seed(u:v,2); %get alternates
          hs = hist(second,1:length(second)); %alternate degree connectivity
          nzhs = find(hs~=0); %find non-zero frequencies
          second_links(gen,nzhs) = hist(hs(nzhs),1:length(nzhs)); %frequency of alternate degree
        end
%         end
      end
      sml = cat_row(sml,sum(mate_links,1)./sum(population));
      sal = cat_row(sal,sum(second_links,1)./sum(population));
%       [nzml iml] = find(mate_links~=0);
%       [nzsl isl] = find(second_links~=0);
%       sml = cat_row(sml,sum(mate_links(nzml,iml),1));
%       sal = cat_row(sal,sum(second_links(nzsl,isl),1));
    end
    end
  end
  %% plot connectedness
  if do_connectedness_plots, 
    figure(8161);
    plot(1:length(sml),mean(sml),'bx'); hold on;
    plot(1:length(sal),mean(sal),'rx'); 
    title(make_title_name(base_name,[]));
    xlabel('degree connectivity');
    ylabel('frequency of degree');
    legend('mate','alternate');
  end
end
end
end

%% Old Code:
%       if do_connectedness_plots==1
%         [nzml iml] = find(mate_links~=0);
%         [nzsl isl] = find(second_links~=0);

%         [tn] = make_title_name(base_name,run_name);
% 
%         figure(8159); 
%         subplot(2,1,1); plot(iml,sum(mate_links(nzml,iml),1),'x');
%         xlabel('number of mate links');  ylabel('frequency of occurence'); title(tn);
% 
%         subplot(2,1,2); plot(isl,sum(second_links(nzsl,isl),1),'x');
%         xlabel('number of second neighbor links');  ylabel('frequency of occurence'); title(tn);
        
%         subplot(2,1,1); plot(iml,log10(sum(mate_links(nzml,iml),1)/(sum(population)*(limit-1))),'x');
%         xlabel('number of mate links');  ylabel('frequency of occurence'); title(tn);
% 
%         subplot(2,1,2); plot(isl,log10(sum(second_links(nzsl,isl),1)/(sum(population)*(limit-1))),'x');
%         xlabel('number of second neighbor links');  ylabel('frequency of occurence'); title(tn);

  %       figure(8160); 
  %       subplot(2,1,1); loglog(iml,sum(mate_links(nzml,iml),1),'x');
  %       xlabel('number of mate links');  ylabel('frequency of occurence'); title(tn);
  % 
  %       subplot(2,1,2); loglog(isl,sum(second_links(nzsl,isl),1),'x');
  %       xlabel('number of second neighbor links');  ylabel('frequency of occurence'); title(tn);
  % 
  %       figure(8161); 
  %       subplot(2,1,1); plot(iml,sum(mate_links(nzml,iml),1),'x');
  %       xlabel('number of mate links');  ylabel('frequency of occurence'); title(tn);
  % 
  %       subplot(2,1,2); plot(isl,sum(second_links(nzsl,isl),1),'x');
  %       xlabel('number of second neighbor links');  ylabel('frequency of occurence'); title(tn);
%       end
global SIMOPTS;
b = 0;
for op = overpop, SIMOPTS.op = op;
for dm = death_max, SIMOPTS.dm = dm;
for mu = mutability, SIMOPTS.mu = mu;
  make_dir = 0; [base_name,dir_name] = NameAndCD(make_dir);
  b = b +1;
  in_degree_frequencies = [];
  min_deg_old = 10^10;
  max_deg_old = 0;
  for run = SIMS
    run_name = int2str(run);
%     if exist([make_data_name('in_degree',base_name,run_name,exp_type,reproduction,limit,0) '.mat'])~=2 
      load(make_data_name('population',base_name,run_name,0));
      load(make_data_name('num_clusters',base_name,run_name,0));
      load(make_data_name('trace_cluster_seed',base_name,run_name,0));
      load(make_data_name('trace_cluster',base_name,run_name,0));
      ngen = length(find(population>=limit));
      mate_links = zeros(ngen,10000,'uint16');
      if size(trace_cluster_seed,2)>=2, second_links = zeros(ngen,10000,'uint16'); 
      else second_links = []; end
      u = 0;  v = 0;
      in_degree = zeros(sum(population),1);
      cluster_coefficient = zeros(sum(population),1);
      for gen = 1:ngen
        u = v +1;
        v = sum(population(1:gen));
        for tci = 1:num_clusters(gen)
          these_guys = find(trace_cluster(u:v)==tci);
          mates = trace_cluster_seed(u-1+these_guys,1);
          if limit>2, alternates = trace_cluster_seed(u-1+these_guys,2);  else, alternates = [];  end
          for tgi = 1:length(these_guys)
            i = u +these_guys(tgi) -1;
            mates_of_tgi = find(mates==these_guys(tgi));
            alternates_of_tgi = find(alternates==these_guys(tgi));
            in_degree(i) = length(mates_of_tgi);
            if limit>2, in_degree(i) = in_degree(i) +length(alternates_of_tgi); end
            total_degree = in_degree(i) +limit -1;
            unique_degree = length(unique([these_guys(mates_of_tgi)', ...
                                           these_guys(alternates_of_tgi)', ...
                                           trace_cluster_seed(tgi,:)]));
            cluster_coefficient(i) = ((limit -1)*length(these_guys))/(total_degree*(total_degree -1)/2);
          end
        end
      end
      in_degree_frequency = hist(in_degree,length(min(in_degree):max(in_degree)));
      min_deg_new = min(in_degree);
      if min_deg_new<min_deg_old
        min_deg = min_deg_new;
        min_deg_old = min_deg_new;
      end
      max_deg_new = max(in_degree);
      if max_deg_new>max_deg_old
        max_deg = max_deg_new;
        max_deg_old = max_deg_new;
      end
      in_degree_frequencies = cat_row(in_degree_frequencies,in_degree_frequency);
      in_deg_name = make_data_name('in_degree',base_name,run_name,0);
      in_deg_freq_name = make_data_name('in_degree_frequency',base_name,run_name,0);
      cc_name = make_data_name('cluster_coefficient',base_name,run_name,0);
      if exist([in_deg_name '.mat'])~=2 %#ok<*EXIST>
        save(in_deg_name,'in_degree');
      end
      if exist([in_deg_freq_name '.mat'])~=2
        save(in_deg_freq_name,'in_degree_frequency');
      end
      if exist([cc_name '.mat'])~=2
        save(cc_name,'cluster_coefficient');
      end
      prev = 0;
%     else
%       load(make_data_name('in_degree_frequency',base_name,run_name,exp_type,reproduction,limit,0));
%       in_degree_frequencies = cat_row(in_degree_frequencies,in_degree_frequency);
%       prev = 1;
%     end
  end
  if do_degrees_plot==1
    [tn] = [make_title_name(base_name,'') 'probability distribution'];
    figure(1510000 +mu*100); 
    midf = mean(in_degree_frequencies,1);
%     stdidf = std(in_degree_frequencies,1);
    if prev==0
      plot(min_deg:max_deg,midf/sum(midf),'x');
    else
      plot(0:length(midf)-1,midf/sum(midf),'x');
    end
    xlabel('in\_degree');  ylabel('in\_degree probability'); title(tn);
  end
end
end
end
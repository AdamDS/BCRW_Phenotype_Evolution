%       [k1] = print_check_result(str_k,ck(1,1),equals,checks_k(1),kerror,ck(2,1));
%       [k2] = print_check_result(str_k,ck(1,2),lte,checks_k(2),kerror,ck(2,2));
%       [k3] = print_check_result(str_k,ck(1,3),gte,checks_k(3),kerror,ck(2,3));
%       kout = [k1,k2,k3];
% %         if gor==1,  
% %           len_r = size(rivalries,1);  
% %           if sum_p==len_r,  c_lenr = 1; end 
% %           min_r = min(rivalries);
% %           if min_k<=min_r,  c_minr = 1; end
% %           max_r = max(rivalries);
%           max_kop = max(kills(:,1));
% %           if max_kop>=max_r,  c_maxr = 1; end
% %         end
% %       end
%       checks_r = [sum_p,min_k,max_kop];
%       [cr] = data_check_B(gor,rivalries,checks_r);
%       [r1] = print_check_result(str_r,cr(1,1),equals,checks_r(1),rerror,cr(2,1));
%       [r2] = print_check_result(str_r,cr(1,2),lte,checks_r(2),rerror,cr(2,2));
%       [r3] = print_check_result(str_r,cr(1,3),gte,checks_r(3),rerror,cr(2,3));
%       rout = [r1,r2,r3];
%       spmi = sum_p -ipop;
% %       if gopar==1,  
% %         len_par = size(parents,1);  
% %         if sum_p==len_par,  c_lenpar = 1; end 
% %         min_par = min(min(parents));
% %         if 1<=min_par,  c_minpar = 1; end
% %         max_par = max(max(parents));
% %         if max_p>=max_p,  c_maxpar = 1; end
% %       end
%       checks_par = [spmi,1,max_p];
%       [cpar] = data_check_B(gopar,parents,checks_par);
%       [par1] = print_check_result(str_par,cpar(1,1),equals,checks_par(1),parerror,cpar(2,1));
%       [par2] = print_check_result(str_par,cpar(1,2),lte,checks_par(2),parerror,cpar(2,2));
%       [par3] = print_check_result(str_par,cpar(1,3),gte,checks_par(3),parerror,cpar(2,3));
%       parout = [par1,par2,par3];
%       %% Check num_clusters dependent dimensions
% %       if gonc==1, 
% %         len_nc = length(num_clusters);  
% %         if len_p==len_nc, c_lennc = 1;  end
% %         min_nc = min(num_clusters);
%         min_poss_nc = floor(min_p/limit);
% %         if min_poss_nc<=min_nc, c_minnc = 1;  end
%         max_nc = max(num_clusters);
%         max_poss_nc = floor(max_p/limit);
% %         if max_poss_nc>=max_nc, c_maxnc = 1;  end
%         checks_nc = [len_p,min_poss_nc,max_poss_nc];
%         [cnc] = data_check_B(gonc,num_clusters,checks_nc);
%         [nc1] = print_check_result(str_nc,cnc(1,1),equals,checks_nc(1),ncerror,cnc(2,1));
%         [nc2] = print_check_result(str_nc,cnc(1,2),lte,checks_nc(2),ncerror,cnc(2,2));
%         [nc3] = print_check_result(str_nc,cnc(1,3),gte,checks_nc(3),ncerror,cnc(2,3));
%         ncout = [nc1,nc2,nc3];

%         checks_tc = [sum_p,1,max_nc];
%         [ctc] = data_check_B(go

%         sum_nc = sum(num_clusters);
%         if gotc==1, 
%           len_tc = length(trace_cluster); 
%           if sum_p==len_tc,  c_lentc = 1;  end 
%           min_tc = min(trace_cluster);
%           if 1==min_tc, c_mintc = 1;  end
%           max_tc = max(trace_cluster);
%           if max_nc==max_tc,  c_maxtc = 1;  end
%         end
%         if goonc==1,  
%           len_onc = length(orgsncluster); 
%           if sum_nc==len_onc, c_lenonc = 1; end 
%           min_onc = min(orgsncluster);
%           if limit<=min_onc,  c_minonc = 1; end
%           max_onc = max(orgsncluater);
%           if max_p>max_onc, c_maxonc = 1; end
%         end
%         
%         if gocx==1,
%           len_cx = length(centroid_x);  
%           if sum_nc==len_cx,  c_lencx = 1;  end 
%           min_cx = min(centroid_x);
%           if min_land<=min_cx && min_tx<=min_cx, c_mincx = 1; end
%           max_cx = max(centroid_x);
%           if max_land>=max_cx && max_tx>=max_cx,  c_maxcx = 1;  end
%         end
%         
%         if gocy==1, 
%           len_cy = length(centroid_y);  
%           if sum_nc==len_cy,  c_lency = 1;  end 
%           min_cy = min(centroid_y);
%           if min_land<=min_cy && min_ty<=min_cy,  c_mincy = 1;  end
%           max_cy = max(centroid_x);
%           if max_land>=max_cy && max_ty>=max_cy,  c_maxcy = 1;  end
%         end
%                 
%         if godiv==1,  
%           len_div = length(cluster_diversity);  
%           if sum_nc==len_div, c_lendiv = 1; end 
%           min_div = min(cluster_diversity);
%           if op<=min_div, c_mindiv = 1; end
%           max_div = max(cluster_diversity);
%           if diag_land>=max_div,  c_maxdiv = 1; end
%         end
%         
%         if godia==1,  
%           len_dia = length(cluster_diameters);  
%           if sum_nc==len_dia, c_lendia = 1; end 
%           min_dia = min(cluster_diameters);
%           if op<=min_dia, c_mindia = 1; end 
%         end
%         
%         if gopl==1, 
%           len_pl = length(path_length); 
%           if sum_nc==len_pl,  c_lenpl = 1;  end 
%           min_pl = min(path_length);
%           min_poss_pl = (limit-1)*limit*op;
%           if min_poss_pl<=min_pl, c_minpl = 1;  end
%         end
%         
%         lpmo = len_p -1;
%         if gond==1, 
%           [row_nd,col_nd] = size(num_descendants); 
%           if lpmo==row_nd,  c_rownd = 1;  end
%           if ipop==col_nd,  c_colnd = 1;  end 
%           min_nd = min(min(num_descendants));
%           if 0<=min_nd, c_minnd = 1;  end
%           max_nd = max(max(num_descendants));
%           if max_p>=max_nd, c_maxnd = 1;  end
%         end
%         
%         if gondc==1,  
%           [row_ndc,col_ndc] = size(num_descendant_clusters);  
%           if lpmo==row_ndc, c_rowndc = 1; end 
%           if ipop==col_ndc, c_colndc = 1; end 
%           min_ndc = min(min(num_descendant_clusters));
%           if 0<=min_ndc,  c_minndc = 1; end
%           max_ndc = max(max(num_descendant_clusters));
%           if max_nc>=max_ndc, c_maxndc = 1; end
%           
%           if godc==1, 
%             len_dc = length(descendant_clusters); 
%             sum_ndc = sum(sum(num_descendant_clusters));
%             if sum_ndc==len_dc, c_lendc = 1;  end
%             min_dc = min(descendant_clusters);
%             if 1<=min_dc, c_mindc = 1;  end
%             max_dc = max(max(descendant_clusters));
%             if max_nc==max_dc,  c_maxdc = 1;  end
%           end
%         end
%         
%         sncmf = sum_nc -num_clusters(end);
%         if goncp==1,  
%           len_ncp = length(num_clusters_produced);  
%           if sncmf==len_ncp, c_lenncp = 1; end 
%           min_ncp = min(min(num_clusters_produced));
%           if 0<=min_ncp,  c_minncp = 1; end
%           max_ncp = max(max(num_clusters_produced));
%           if max_nc>=max_ncp, c_maxncp = 1; end
%           
%           sum_ncp = sum(num_clusters_produced);
%           if gocp==1,  
%             len_cp = length(clusters_produced);  
%             if sum_ncp==len_cp, c_lencp = 1; end 
%             min_cp = min(clusters_produced);
%             if 1<=min_cp, c_mincp = 1;  end
%             max_cp = max(clusters_produced);
%             if max_nc==max_cp,  c_maxcp = 1;  end
%           end
%         end
%         
%         sncmi = sum_nc -num_clusters(1);
%         if goncf==1,  
%           len_ncf = length(num_clusters_fused);  if sncmi==len_ncf, c_lenncf = 1; end 
%           sum_ncf = sum(num_clusters_fused);
%           if gocf==1,  
%             len_cf = length(clusters_fused);  if sum_ncf==len_cf, c_lencf = 1; end
%           end
%         end
%       end
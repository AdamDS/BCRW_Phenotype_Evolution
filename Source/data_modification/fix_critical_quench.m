usual = [1:6 11:15 21:25 31:35 42:45 51:55];

long = zeros(100,1);
long(usual) = 1;

NGEN = 10^5;
pname = 'populations_Bacterial_Mu_28_Uniform_2_Flatscape_10x10_basic_map_size_20290_IPOP_100000_NGEN_1_100.mat';
load(pname);

if size(populations,2)>NGEN,  
  figure, plot(mean(populations(long==1,:)),'r')
  hold on;
  plot(mean(populations(long==0,:)),'b')
  plot(mean(populations),'k')
  set(gca,'xscale','log','yscale','log');
  cutoff = input('enter cutoff gen');
  for i=1:size(populations,1), lt(i) = length(find(populations(i,1:cutoff))); end
  mlt = max(lt);
  fprintf('max lifetime %d\n',mlt);
end
yn = input('save fixed populations? y/n');
if yn,  
  save(pname,'populations');
end
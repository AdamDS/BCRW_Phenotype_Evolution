function [] = save_lifetimes(base_name,dir_name)
global lifetimes;
lt_name = [dir_name 'lifetimes_' base_name];
save(lt_name,'lifetimes');
end
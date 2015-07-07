%% data_check_.m
%[couts] = data_check_(flag,data,compval)
%

function [couts] = data_check_(data,compval), 
couts = zeros(2,3);
initial = data(1);
if compval(1)==initial, cout(1) = 1;  end
min_ = min(data);
if min_>=compval(2), cout(2) = 1; end
max_ = max(data);
len_ = length(data);
find_ = find(data>=compval(2));
if len_==length(find_), cout(3) = 1; end
end
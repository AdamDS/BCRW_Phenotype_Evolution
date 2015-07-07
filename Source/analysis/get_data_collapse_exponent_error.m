%% get_data_collapse_exponent_error.m
%
% [dexp] = get_data_collapse_exponent_error(best_exp,exponents,eta,Pb,best_Pb,sgn)
%
function [dexp] = get_data_collapse_exponent_error(best_exp,exponents,eta,Pb,best_Pb,sgn),  
oode = 1/mean(diff(exponents)); %One Over Difference between Exponents
plus = (1+sgn*eta)*best_exp;
rex = abs(exponents-(round(plus*oode)/oode));
Pbp = Pb(rex==min(rex(rex~=0))); %get nearest Pb(exponents)
dexp = eta*best_exp*sqrt(2*log(Pbp(1)/best_Pb)).^-1;
end
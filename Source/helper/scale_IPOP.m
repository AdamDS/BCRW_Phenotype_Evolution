function [IPOP] = scale_IPOP(),  
global SIMOPTS;
maxPop = get_max_population(); 
den = SIMOPTS.defIPOP/37544;
IPOP = floor(den*maxPop);
SIMOPTS.IPOP = IPOP;
end
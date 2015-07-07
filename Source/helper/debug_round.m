%% debug_round.m
% [] = debug_round(counter)
%
% How to use:
% 1) Initialize a counter at the beginning of a function. (counter = 0)
% 2) Place debug_round(counter) before the line where an error is thrown.
% 3) Note the last counter value in the command window.
% 4) Place a breakpoint on the error throwing line, and set a conditional
% for when counter==last counter value from the command window.
% 5) Figure out the bug just before it's thrown and fix it if possible.
% 6) Remove counter initialization and debug_round from code.
% 7) Enjoy!
%
% 8*15*13 ADS

function [counter] = debug_round(counter), 
counter = counter +1;
fprintf(['debug round %1.0f \n'],counter);
end
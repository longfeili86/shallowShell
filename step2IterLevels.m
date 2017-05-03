function [prev2,prev,cur,new] = step2IterLevels(step)
% convert step to the corresponding iteration level.
numberOfLevels=4; % we save four levels of solutions, we need for levels to estimate conv rate
cur=mod(step,numberOfLevels)+1; % matlab index starts from, so add 1
new=mod(step+1,numberOfLevels)+1;
prev=mod(step-1,numberOfLevels)+1;
prev2=mod(step-2,numberOfLevels)+1;

end
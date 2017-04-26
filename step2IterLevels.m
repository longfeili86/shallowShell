function [prev,cur,new] = step2IterLevels(step)
% convert step to the corresponding iteration level.
numberOfLevels=3; % we save three levels of solutions
cur=mod(step,numberOfLevels)+1; % matlab index starts from, so add 1
new=mod(step+1,numberOfLevels)+1;
prev=mod(step-1,numberOfLevels)+1;

end
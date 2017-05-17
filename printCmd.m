function printCmd(options)    
% this function print out the command if the options are specified as a
% cell

fprintf('runShell ');
for i=1:length(options)
    fprintf('%s ',options{i});
end
fprintf('\n');
    
end
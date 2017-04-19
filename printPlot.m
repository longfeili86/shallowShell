function printPlot(figureName,resultsDir)
% this function print the plots to a eps file to resultsDir
% -- Longfei Li

if(nargin==1)
    resultsDir='.'; % default is pwd
end
infoPrefix='---printPlot--: ';    

figureName=sprintf('%s/%s',resultsDir,figureName);
fprintf('%splot saved. filename=%s\n',infoPrefix,figureName);
print('-depsc2',figureName);

end
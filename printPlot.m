function printPlot(resultsDir,figureName)
infoPrefix='---printPlot--: ';    

figureName=sprintf('%s/%s',resultsDir,figureName);
fprintf('%splot saved. filename=%s\n',infoPrefix,figureName);
print('-depsc2',figureName);

end
function moveFigureInSavedResults(resultsDir,destination)
% move eps files in resultsDir to destination
% eps was renamed with resultsDir as its prefix

fprintf('mv figures from %s to %s\n', resultsDir,destination);
files=dir(sprintf('%s/*.eps',resultsDir));
for i = 1: length(files)
    oldName=sprintf('%s/%s',resultsDir,files(i).name);
    newName=sprintf('%s/%s_%s',destination,resultsDir,files(i).name);
    %fprintf('mv %s %s\n', oldName, newName);
    movefile(oldName,newName);
end



end
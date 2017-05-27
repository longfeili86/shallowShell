function x0=readICfromSavedResults(xiNew,resultsDir1,resultsDir2)
infoPrefix = '--readICfromSavedResults--: '; % all info displayed by this function includes this prefix

assert(exist(resultsDir1,'dir')~=0,...
    'Error: the saved results that provide initial guess does not exist. It needs to be generated first');
if(~strcmp(resultsDir2,''))
    assert(exist(resultsDir2,'dir')~=0,...
        'Error: the second saved results that provide initial guess does not exist. It needs to be generated first');

    % use two levels of previous 
    fprintf('%sOtain initial guess using TWO previously saved results: %s and %s\n',infoPrefix,resultsDir1,resultsDir2);
    load(sprintf('%s/results.mat',resultsDir1),'x','xi');
    xOld=x;
    xiOld=xi;

    load(sprintf('%s/results.mat',resultsDir2),'x','xi');
    xOld2=x;
    xiOld2=xi;

    x0= (xOld-xOld2)/(xiOld-xiOld2)*(xiNew-xiOld)+xOld;

else
    % use previous results
    fprintf('%sObtain initial guess using previously saved results: %s\n',infoPrefix,resultsDir1);
    load(sprintf('%s/results.mat',resultsDir1),'x');
    x0=x;
end


end
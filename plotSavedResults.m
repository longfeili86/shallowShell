function plotSavedResults(varargin)
% this function plot the saved results in resultsDir
% usage: plotSavedResults -f=<resultsDir> (optional: -contour -nosave)
% --Longfei Li
resultsDir='.'; % default dir is pwd
savePlot=true;
contour=false;
for i=1:nargin
    line = varargin{i};
    if(strncmp(line,'-f=',3))
        resultsDir=line(4:end); 
    elseif(strcmp(line,'-nosave'))
        savePlot=false;
    elseif(strcmp(line,'-contour'))
        contour=true;
    end
end

fprintf('Plot saved results from the resultsDir=%s/\n',resultsDir);

load(sprintf('%s/results.mat',resultsDir));
suffix='';
if(contour)
    suffix='Contour';
end

figure
mySurf(Xplot,Yplot,Wplot,'Solution',contour);
if(savePlot)
    printPlot(sprintf('solution%s',suffix),resultsDir);
end

if(exist('Fplot','var'))
    figure
    mySurf(Xplot,Yplot,Fplot,'Forcing',contour);
    if(savePlot)
        printPlot(sprintf('forcing%s',suffix),resultsDir);
    end
end

%plot error if exist
if(exist('errPlot','var'))
    figure
    mySurf(Xplot,Yplot,errPlot,'Error',contour);
    if(savePlot)
        printPlot(sprintf('error%s',suffix),resultsDir);
    end    
end



end
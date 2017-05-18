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
mySurf(Xplot,Yplot,Wplot,'$w$',contour);
if(savePlot)
    printPlot(sprintf('wSolution%s',suffix),resultsDir);
end

if(exist('PHIplot','var'))
    figure
    mySurf(Xplot,Yplot,PHIplot,'$\phi$',contour);
    if(savePlot)
        printPlot(sprintf('phiSolution%s',suffix),resultsDir);
    end
end

if(exist('W0plot','var'))
    figure
    mySurf(Xplot,Yplot,W0plot,'$w_0$',contour);
    if(savePlot)
        printPlot(sprintf('w0%s',suffix),resultsDir);
    end
    % plot w+w0
    figure
    mySurf(Xplot,Yplot,Wplot+W0plot,'$w+w_0$',contour);
    if(savePlot)
        printPlot(sprintf('wPlusw0%s',suffix),resultsDir);
    end        
end

if(exist('Fwplot','var'))
    figure
    mySurf(Xplot,Yplot,Fwplot,'$f_w$',contour);
    if(savePlot)
        printPlot(sprintf('wForcing%s',suffix),resultsDir);
    end
end

if(exist('Fphiplot','var'))
    figure
    mySurf(Xplot,Yplot,Fphiplot,'$f_{\phi}$',contour);
    if(savePlot)
        printPlot(sprintf('phiForcing%s',suffix),resultsDir);
    end
end

%plot error if exist
if(exist('WerrPlot','var'))
    figure
    mySurf(Xplot,Yplot,WerrPlot,'$E(w)$',contour);
    if(savePlot)
        printPlot(sprintf('wError%s',suffix),resultsDir);
    end    
end

if(exist('PHIerrPlot','var'))
    figure
    mySurf(Xplot,Yplot,PHIerrPlot,'$E(\phi)$',contour);
    if(savePlot)
        printPlot(sprintf('phiError%s',suffix),resultsDir);
    end    
end


end
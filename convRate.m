function convRate(varargin)
%this function does convergence study
%--------------- some default options -----------------
% currently available tests: 'biharmPolyTest','biharmTrigTest'; 
% Usage: convRate -case=<caseName> -test=<testName> -rStart=<int> -rEnd=<int> (optional -run -noplot)
% -- Longfei Li

%print current time
fprintf('%s\n',datestr(now));

% default parameters
testName='biharmTrigTest'; % available tests are now: biharmTrigTest biharmPolyTest
caseName='biharmonic'; % do conv study for biharmonic solver
solver=''; % solver for coupledSystem: fsolve, exPicard, imPicard
isRun=false; % flag to do the computation, otherwise just do the plots from saved data
isPlot=true;
rStart=1; % refinement start number
rEnd=4; % refinement end number

for i=1:nargin
    line = varargin{i};
    if(strncmp(line,'-test=',6))
        testName=line(7:end); 
    elseif(strncmp(line,'-case=',6))
        caseName=line(7:end);
    elseif(strncmp(line,'-solver=',8))
        solver=line(9:end); 
    elseif(strncmp(line,'-rStart=',8))
        rStart=sscanf(line,'-rStart=%i'); 
    elseif(strncmp(line,'-rEnd=',6))
        rEnd=sscanf(line,'-rEnd=%i');        
    elseif(strcmp(line,'-run'))
        isRun=true;  
    elseif(strcmp(line,'-noplot'))
        isPlot=false;         
    end
end


rf=rStart:rEnd; % do conv study on these refinements
bcNames={'Supported','Clamped','Free'};% avaible bcs


gn=@(i) 2.^(rf(i));
numGrids=@(i) gn(i).*10;
if(strcmp(caseName,'coupledSystem'))
    resultsDir=@(i,bcType,testName) sprintf('%s%s%sG%i',solver,testName,bcNames{bcType},gn(i));
else
    resultsDir=@(i,bcType,testName) sprintf('%s%sG%i',testName,bcNames{bcType},gn(i));
end

% call runShell to do the computation on a serie of refined grids
if(isRun)    
    funcDefFilename=sprintf('%sFuncDef',testName);
  
    % define cmd 
    if(strcmp(caseName,'coupledSystem'))
        getCMD = @(i,bcType,testName) sprintf('runShell -case=%s -solver=%s -funcDefFile=%s -nx=%i -ny=%i -bcType=%i -noplot -f=%s',...
            caseName,solver,funcDefFilename,numGrids(i),numGrids(i),bcType,resultsDir(i,bcType,testName));        
    else
        getCMD = @(i,bcType,testName) sprintf('runShell -case=%s -funcDefFile=%s -nx=%i -ny=%i -bcType=%i -noplot -f=%s',...
            caseName,funcDefFilename,numGrids(i),numGrids(i),bcType,resultsDir(i,bcType,testName));
    end
    for bcType=1:3
        for i=1:length(rf)
            cmd=getCMD(i,bcType,testName);
            fprintf('--running--: %s\n',cmd);
            eval(cmd);
        end
    end
end

% plot convRate
if(isPlot)    
    % read data 
    for bcType=1:3
        for i=1:length(rf)
            load(sprintf('%s/results.mat',resultsDir(i,bcType,testName)));
            wErrorL2.(bcNames{bcType})(i)=norm(WerrPlot(:),2)/sqrt(length(WerrPlot(:)));
            wErrorLmax.(bcNames{bcType})(i)=max(abs(WerrPlot(:)));
            if(strcmp(caseName,'coupledSystem'))
                phiErrorL2.(bcNames{bcType})(i)=norm(PHIerrPlot(:),2)/sqrt(length(PHIerrPlot(:)));
                phiErrorLmax.(bcNames{bcType})(i)=max(abs(PHIerrPlot(:)));    
            end
        end
    end  

    hh=1./numGrids(1:length(rf));
    colors={'r','g','b','c','m'};
    
    % Lmax convRate
    figure
    setupFigure;
    loglog(hh,(hh).^2,'k','LineWidth',figOptions.LW,'MarkerSize',figOptions.MS);
    hold on
        for bcType=1:3
        loglog(hh,wErrorLmax.(bcNames{bcType}),strcat(colors{bcType},'-+'),'LineWidth',figOptions.LW,'MarkerSize',figOptions.MS);
        lgdNames{bcType}=strcat(bcNames{bcType},' ($w$)');
        end
        if(strcmp(caseName,'coupledSystem'))
            for bcType=1:3
            loglog(hh,phiErrorLmax.(bcNames{bcType}),strcat(colors{bcType},'->'),'LineWidth',figOptions.LW,'MarkerSize',figOptions.MS); 
            lgdNames{bcType+3}=strcat(bcNames{bcType},' ($\phi$)');
            end
        end
    hold off
    legend(['2nd order',lgdNames],'Location','NorthWest','FontSize',figOptions.FS,'Interpreter','Latex');
    xlabel('grid size','FontSize',figOptions.FS) ;
    ylabel('error','FontSize',figOptions.FS) ;
    set(gca,'FontSize',figOptions.FS);
    box on
    grid on
    if(strcmp(caseName,'coupledSystem'))
        figName=sprintf('%s%sConvRateLmax.eps',solver,testName);
    else
        figName=sprintf('%sConvRateLmax.eps',testName);
    end
    printPlot(figName);
    
    % L2 convRate
    figure
    setupFigure;
    loglog(hh,(hh).^2,'k','LineWidth',figOptions.LW,'MarkerSize',figOptions.MS);
    hold on
        for bcType=1:3
        loglog(hh,wErrorL2.(bcNames{bcType}),strcat(colors{bcType},'-+'),'LineWidth',figOptions.LW,'MarkerSize',figOptions.MS);
        lgdNames{bcType}=strcat(bcNames{bcType},' ($w$)');
        end
        if(strcmp(caseName,'coupledSystem'))
            for bcType=1:3
            loglog(hh,phiErrorL2.(bcNames{bcType}),strcat(colors{bcType},'->'),'LineWidth',figOptions.LW,'MarkerSize',figOptions.MS); 
            lgdNames{bcType+3}=strcat(bcNames{bcType},' ($\phi$)');
            end
        end
    hold off
    legend(['2nd order',lgdNames],'Location','NorthWest','FontSize',figOptions.FS,'Interpreter','Latex');
    xlabel('grid size','FontSize',figOptions.FS) ;
    ylabel('error','FontSize',figOptions.FS) ;
    set(gca,'FontSize',figOptions.FS);
    box on
    grid on
    if(strcmp(caseName,'coupledSystem'))
        figName=sprintf('%s%sConvRateL2.eps',solver,testName);
    else
        figName=sprintf('%sConvRateL2.eps',testName);
    end
    printPlot(figName);
end

end




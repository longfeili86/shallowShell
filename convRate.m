function convRate(varargin)
%this function does convergence study
%--------------- some default options -----------------
% currently available tests: 'biharmPolyTest','biharmTrigTest'; 
% Usage: convRate -case=<caseName> -test=<testName>  (optional -run -noplot)
% -- Longfei Li

%print current time
fprintf('%s\n',datestr(now));

% default parameters
testName='biharmTrigTest'; % available tests are now: biharmTrigTest biharmPolyTest
caseName='biharmonic'; % do conv study for biharmonic solver
isRun=false; % flag to do the computation, otherwise just do the plots from saved data
isPlot=true;

for i=1:nargin
    line = varargin{i};
    if(strncmp(line,'-test=',6))
        testName=line(7:end); 
    elseif(strncmp(line,'-case=',6))
        caseName=line(7:end);        
    elseif(strcmp(line,'-run'))
        isRun=true;  
    elseif(strcmp(line,'-noplot'))
        isPlot=false;         
    end
end


rf=2:7; % do conv study on these refinements
bcNames={'Supported','Clamped','Free'};% avaible bcs


gn=@(i) 2.^(rf(i)-1);
numGrids=@(i) gn(i).*10;
resultsDir=@(i,bcType,testName) sprintf('%s%sG%i',testName,bcNames{bcType},gn(i));

% call runShell to do the computation on a serie of refined grids
if(isRun)    
    rhsFilename=sprintf('%sRHS',testName);
  
    % define cmd 
    getCMD = @(i,bcType,testName) sprintf('runShell -case=%s -rhsFile=%s -nx=%i -ny=%i -bcType=%i -noplot -f=%s',...
        caseName,rhsFilename,numGrids(i),numGrids(i),bcType,resultsDir(i,bcType,testName));
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
            ErrorL2.(bcNames{bcType})(i)=norm(errPlot(:),2)/sqrt(length(errPlot(:)));
            ErrorLmax.(bcNames{bcType})(i)=max(abs(errPlot(:)));
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
        loglog(hh,ErrorLmax.(bcNames{bcType}),strcat(colors{bcType},'-+'),'LineWidth',figOptions.LW,'MarkerSize',figOptions.MS); 
        end
    hold off
    legend(['2nd order',bcNames],'Location','NorthWest','FontSize',figOptions.FS);
    xlabel('grid size','FontSize',figOptions.FS) ;
    ylabel('error','FontSize',figOptions.FS) ;
    set(gca,'FontSize',figOptions.FS);
    box on
    grid on
    printPlot(sprintf('%sConvRateLmax.eps',testName));
    
    % L2 convRate
    figure
    setupFigure;
    loglog(hh,(hh).^2,'k','LineWidth',figOptions.LW,'MarkerSize',figOptions.MS);
    hold on
        for bcType=1:3
        loglog(hh,ErrorL2.(bcNames{bcType}),strcat(colors{bcType},'-+'),'LineWidth',figOptions.LW,'MarkerSize',figOptions.MS); 
        end
    hold off
    legend(['2nd order',bcNames],'Location','NorthWest','FontSize',figOptions.FS);
    xlabel('grid size','FontSize',figOptions.FS) ;
    ylabel('error','FontSize',figOptions.FS) ;
    set(gca,'FontSize',figOptions.FS);
    box on
    grid on
    printPlot(sprintf('%sConvRateL2.eps',testName));
end

end




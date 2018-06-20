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
isLinear=true;
rStart=1; % refinement start number
rEnd=4; % refinement end number
legendPosition = 'southeast' ; % legend position

for i=1:nargin
    line = varargin{i};
    if(strncmp(line,'-test=',6))
        testName=line(7:end); 
    elseif(strncmp(line,'-case=',6))
        caseName=line(7:end);
    elseif(strncmp(line,'-solver=',8))
        solver=line(9:end);         
    elseif(strcmp(line,'-run'))
        isRun=true;  
    elseif(strcmp(line,'-noplot'))
        isPlot=false; 
    elseif(strcmp(line,'-nonlinear'))
        isLinear=false;        
    elseif(strncmp(line,'-rStart=',8))
        rStart=sscanf(line,'-rStart=%i'); 
    elseif(strncmp(line,'-rEnd=',6))
        rEnd=sscanf(line,'-rEnd=%i');       
    end
end


rf=rStart:rEnd; % do conv study on these refinements
bcNames={'Supported','Clamped','Free','CS','CF'};% avaible bcs
numBCTypes=5;
bcTypes=1:numBCTypes; % do conv test for these bcTypes


gn=@(i) 2.^(rf(i));
numGrids=@(i) gn(i).*10;
if(strcmp(caseName,'coupledSystem'))
    linearity='Nonlinear';
    if(isLinear)
        linearity='Linear';
    end
    resultsDir=@(i,bcType,testName) sprintf('%s%s%s%sG%i',solver,linearity,testName,bcNames{bcType},gn(i));
    %bcTypes=setdiff(bcTypes,3); % free bc needs more resolution. We deal it separately
else
    resultsDir=@(i,bcType,testName) sprintf('%s%sG%i',testName,bcNames{bcType},gn(i));
end

% call runShell to do the computation on a serie of refined grids
if(isRun)    
    funcDefFilename=sprintf('%sFuncDef',testName);
  
    % define cmd 
    if(strcmp(caseName,'coupledSystem'))
        if(isLinear)
            getCMD = @(i,bcType,testName) sprintf('runShell -case=%s -solver=%s -funcDefFile=%s -nx=%i -ny=%i -bcType=%i -noplot -f=%s',...
                caseName,solver,funcDefFilename,numGrids(i),numGrids(i),bcType,resultsDir(i,bcType,testName)); 
        else
            getCMD = @(i,bcType,testName) sprintf('runShell -case=%s -solver=%s -funcDefFile=%s -nx=%i -ny=%i -bcType=%i -noplot -f=%s -nonlinear',...
                caseName,solver,funcDefFilename,numGrids(i),numGrids(i),bcType,resultsDir(i,bcType,testName));            
        end
    else
        getCMD = @(i,bcType,testName) sprintf('runShell -case=%s -funcDefFile=%s -nx=%i -ny=%i -bcType=%i -noplot -f=%s',...
            caseName,funcDefFilename,numGrids(i),numGrids(i),bcType,resultsDir(i,bcType,testName));
    end
    for bcType=bcTypes
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
    for bcType=bcTypes
        for i=1:length(rf)
            load(sprintf('%s/results.mat',resultsDir(i,bcType,testName)));
            wErrorL2.(bcNames{bcType})(i)=norm(WerrPlot(:),2)/sqrt(length(WerrPlot(:)));
            wErrorLmax.(bcNames{bcType})(i)=max(abs(WerrPlot(:)));
                
            if(strcmp(caseName,'coupledSystem'))
                phiErrorL2.(bcNames{bcType})(i)=norm(PHIerrPlot(:),2)/sqrt(length(PHIerrPlot(:)));
                phiErrorLmax.(bcNames{bcType})(i)=max(abs(PHIerrPlot(:)));  
                if(wErrorLmax.(bcNames{bcType})(i)>1e1 || phiErrorLmax.(bcNames{bcType})(i)>1e1)
                     % chop off unconverged results
                    wErrorL2.(bcNames{bcType})(i)=nan;
                    wErrorLmax.(bcNames{bcType})(i)=nan;    
                    phiErrorL2.(bcNames{bcType})(i)=nan;
                    phiErrorLmax.(bcNames{bcType})(i)=nan;
                end
            end
        end
    end  

    hh=1./numGrids(1:length(rf));
    colors={'b','r','g','c','m'};
    
    % Lmax convRate
    figure
    setupFigure;
    loglog(hh,(hh).^2*3,'k','LineWidth',figOptions.LW,'MarkerSize',figOptions.MS);
    hold on
        counter=0;
        for bcType=bcTypes
        loglog(hh,wErrorLmax.(bcNames{bcType}),strcat(colors{bcType},'-+'),'LineWidth',figOptions.LW,'MarkerSize',figOptions.MS);
        counter=counter+1;
        lgdNames{counter}=strcat(bcNames{bcType},' ($w$)');
        end
        if(strcmp(caseName,'coupledSystem'))
            for bcType=bcTypes
            loglog(hh,phiErrorLmax.(bcNames{bcType}),strcat(colors{bcType},'->'),'LineWidth',figOptions.LW,'MarkerSize',figOptions.MS);
            counter=counter+1;
            lgdNames{counter}=strcat(bcNames{bcType},' ($\phi$)');
            end
        end
    hold off
    legend(['2nd order',lgdNames],'Location',legendPosition,'FontSize',figOptions.FS,'Interpreter','Latex');
    xlabel('grid size','FontSize',figOptions.FS) ;
    ylabel('error','FontSize',figOptions.FS) ;
    set(gca,'FontSize',figOptions.FS);
    box on
    grid off
    xlim([1e-3,10^(-0.2)]);
    ylim([1e-6,1e1]);
    if(strcmp(caseName,'coupledSystem'))
        figName=sprintf('%s%s%sConvRateLmax.eps',solver,linearity,testName);
    else
        figName=sprintf('%sConvRateLmax.eps',testName);
    end
    printPlot(figName);
    
    % L2 convRate
    figure
    setupFigure;
    loglog(hh,(hh).^2*3,'k','LineWidth',figOptions.LW,'MarkerSize',figOptions.MS);
    hold on
        counter=0;
        for bcType=bcTypes
        loglog(hh,wErrorL2.(bcNames{bcType}),strcat(colors{bcType},'-+'),'LineWidth',figOptions.LW,'MarkerSize',figOptions.MS);
        counter=counter+1;
        lgdNames{counter}=strcat(bcNames{bcType},' ($w$)');
        end
        if(strcmp(caseName,'coupledSystem'))
            for bcType=bcTypes
            loglog(hh,phiErrorL2.(bcNames{bcType}),strcat(colors{bcType},'->'),'LineWidth',figOptions.LW,'MarkerSize',figOptions.MS); 
            counter=counter+1;
            lgdNames{counter}=strcat(bcNames{bcType},' ($\phi$)');
            end
        end
    hold off
    legend(['2nd order',lgdNames],'Location',legendPosition,'FontSize',figOptions.FS,'Interpreter','Latex');
    xlabel('grid size','FontSize',figOptions.FS) ;
    ylabel('error','FontSize',figOptions.FS) ;
    set(gca,'FontSize',figOptions.FS);
    box on
    grid off
    xlim([1e-3,10^(-0.2)]);
    ylim([1e-6,1e1]);
    if(strcmp(caseName,'coupledSystem'))
        figName=sprintf('%s%s%sConvRateL2.eps',solver,linearity,testName);
    else
        figName=sprintf('%sConvRateL2.eps',testName);
    end
    printPlot(figName);
end

end




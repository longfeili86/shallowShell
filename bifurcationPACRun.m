function bifurcationPACRun(varargin)
% this function do the bifurcation computation for various cases
infoPrefix='--bifurcationRun-- ';
% default parameters
bcType=1;
isPlot=true;
nx=40;
ny=40;
xiMin=-1650.;
xiMax=0.;
dxi=100;
maxIter=100;
tol=1e-6;
relaxFactor1=1.; % relaxation only implemented for newton. This value won't be seen by other methods
relaxFactor2=1.;
relaxFactor3=1.;
implicitFactor1=1.; % implicitFactor only for imPicard method. 1 for fully implicit, 0 for explicit
implicitFactor2=1.;
implicitFactor3=1.;
solver1='newton'; %solver for top branch
solver2='newton'; %solver for middle branch
solver3='newton'; %solver for lower branch
resultsName='bifurcationPACTest';
funcDefFile='bifurcationFuncDef';
for i=1:nargin
    line = varargin{i};
    if(strncmp(line,'-bcType=',7))
        bcType=sscanf(line,'-bcType=%i');
    elseif(strcmp(line,'-noplot'))
        isPlot=false;
    elseif(strncmp(line,'-nx=',4))
        nx=sscanf(line,'-nx=%i');
    elseif(strncmp(line,'-ny=',4))
        ny=sscanf(line,'-ny=%i');
    elseif(strncmp(line,'-xiMin=',7))
        xiMin=sscanf(line,'-xiMin=%e'); 
    elseif(strncmp(line,'-xiMax=',7))
        xiMax=sscanf(line,'-xiMax=%e');  
    elseif(strncmp(line,'-dxi=',5))
        dxi=sscanf(line,'-dxi=%e'); 
    elseif(strncmp(line,'-maxIter=',9))
        maxIter=sscanf(line,'-maxIter=%i');
    elseif(strncmp(line,'-tol=',5))
        tol=sscanf(line,'-tol=%e'); 
    elseif(strncmp(line,'-relaxFactor1=',13))
        relaxFactor1=sscanf(line,'-relaxFactor1=%e');
    elseif(strncmp(line,'-relaxFactor2=',13))
        relaxFactor2=sscanf(line,'-relaxFactor2=%e');
    elseif(strncmp(line,'-relaxFactor3=',13))
        relaxFactor3=sscanf(line,'-relaxFactor3=%e');     
    elseif(strncmp(line,'-implicitFactor1=',16))
        implicitFactor1=sscanf(line,'-implicitFactor1=%e'); 
    elseif(strncmp(line,'-implicitFactor2=',16))
        implicitFactor2=sscanf(line,'-implicitFactor2=%e');  
    elseif(strncmp(line,'-implicitFactor3=',16))
        implicitFactor3=sscanf(line,'-implicitFactor3=%e');    
    elseif(strncmp(line,'-solver1=',9))
        solver1=line(10:end);
    elseif(strncmp(line,'-solver2=',9))
        solver2=line(10:end);
    elseif(strncmp(line,'-solver3=',9))
        solver3=line(10:end);
    elseif(strncmp(line,'-results=',9))
        resultsName=line(10:end);
    elseif(strncmp(line,'-funcDefFile=',13))
        funcDefFile=sscanf(line,'-funcDefFile=%s');    
    end
end

branches=['a','b','c'];

getOptions=@(solver,implicitFactor,relaxFactor,tol,bcType,nx,ny,maxIter,resultsName,counter,branch,funcDefFile)...
        {'-case=coupledSystem',...
         '-nonlinear',...
         '-saveIC',...
         '-noplot',...
         sprintf('-solver=%s',solver),...
         sprintf('-implicitFactor=%e',implicitFactor),...
         sprintf('-relaxFactor=%e',relaxFactor),...
         sprintf('-tol=%e',tol),...
         sprintf('-bcType=%i',bcType),...
         sprintf('-nx=%i',nx),...
         sprintf('-ny=%i',ny),...
         sprintf('-maxIter=%i',maxIter),...
         sprintf('-f=%s%i%c',resultsName,counter,branch),...
         sprintf('-funcDefFile=%s',funcDefFile),...
         };


% 1st branch (top): from xiMax -> xiMin with step -dxi
counter=0;
ds=dxi;
%ThermalLoading=xiMax-dxi*counter;
exitflag=1;
b=1;
xi=xiMax;
while(xi>=xiMin && exitflag>0)
    counter=counter+1;
    options=getOptions(solver1,implicitFactor1,relaxFactor1,tol,bcType,nx,ny,maxIter,resultsName,counter,branches(b),funcDefFile);

    if(counter>2) % use two previous results as initial guess
        options=[options,...
            {sprintf('-readICResult1=%s',sprintf('%s%i%c',resultsName,counter-1,branches(b))),...
             sprintf('-readICResult2=%s',sprintf('%s%i%c',resultsName,counter-2,branches(b)))},...
             sprintf('-ds=%e',ds),...
             '-usePAC'...
            ];
    else
        ThermalLoading=xiMax-dxi*(counter-1);   
        options=[options,sprintf('-xi=%e',ThermalLoading)];
    end
    printCmd(options);
    exitflag = runShell(options{:});
    if(exitflag<=0)
        % this run is not valid. Do not count it and rm the results dir
        resultDir=sprintf('%s%i%c',resultsName,counter,branches(b));
        fprintf('%sThis run is not valid. Do not count it and rm the results dir: %s\n',infoPrefix,resultDir)
        rmdir(resultDir,'s');
        counter=counter-1;
        continue;
    end
    load(sprintf('%s%i%c/results.mat',resultsName,counter,branches(b)),'xi');
end

% 2nd branch: from xiMin -> xiMax with step dxi
counter=0;
ds=dxi;
exitflag=1;
b=2;
xi=xiMin;
reduceTimes=0; % counting reduce times.
reduceMax=10; % max number of time allowed to reduce ds;
counterReduced=-999; % counter for the number of computations using reduced ds. if negative, means ds is not reduced
reduceFactor=10;
while(xi>=xiMin && exitflag>0 && reduceTimes<reduceMax)
    counter=counter+1;
    options=getOptions(solver1,implicitFactor1,relaxFactor1,tol,bcType,nx,ny,maxIter,resultsName,counter,branches(b),funcDefFile);

    if(counter>2) % use two previous results as initial guess
        options=[options,...
            {sprintf('-readICResult1=%s',sprintf('%s%i%c',resultsName,counter-1,branches(b))),...
             sprintf('-readICResult2=%s',sprintf('%s%i%c',resultsName,counter-2,branches(b)))},...
             sprintf('-ds=%e',ds),...
             '-usePAC'...
            ];
    else
        ThermalLoading=xiMin+dxi*(counter-1);   
        options=[options,sprintf('-xi=%e',ThermalLoading)];
    end
    printCmd(options);
    exitflag = runShell(options{:});
    if(exitflag<=0)
        % this run is not valid. Do not count it and rm the results dir
        resultDir=sprintf('%s%i%c',resultsName,counter,branches(b));
        fprintf('%sThis run is not valid. Do not count it and rm the results dir: %s\n',infoPrefix,resultDir)
        rmdir(resultDir,'s');
        counter=counter-1; 
        if(exitflag==0) % maxIter reached try reduce ds by half
            ds=ds/reduceFactor;
            fprintf('$^&*&^$^&* reducing ds to %e\n',ds);
            pause(0.1);
            exitflag=1;  % change exitflag so we can continue the iteration with smaller ds
            reduceTimes=reduceTimes+1;
            counterReduced=0; % start counting the number of computatations using reduced ds
        end
        continue;
    end
    if(counterReduced>=0)
        counterReduced=counterReduced+1;
    end
    if(counterReduced>5)
        reduceTimes=reduceTimes-1;
        ds=ds*reduceFactor; % gradually going back to original ds
        if(ds==dxi)
            counterReduced=-999; % ds is now back to original
            reduceTimes=0;
        end
    end
    load(sprintf('%s%i%c/results.mat',resultsName,counter,branches(b)),'xi');
end


if(isPlot)
    plotBifurcation(strcat('-f=',resultsName),'-color=b','-distBranch','-showIC');
end

end



























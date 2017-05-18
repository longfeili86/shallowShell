function bifurcationRun(varargin)
% this function do the bifurcation computation for various cases
infoPrefix='--bifurcationRun-- ';
% default parameters
bcType=1;
nx=40;
ny=40;
xiMin=-1650.;
xiMax=0.;
dxi=100;
maxIter=100;
tol=1e-6;
relaxFactor=1.; % relaxation only implemented for newton. This value won't be seen by other methods
solver1='fsolve'; %solver for top branch
solver2='fsolve'; %solver for middle branch
solver3='fsolve'; %solver for lower branch
resultsName='bifurcationTest';
funcDefFile='bifurcationFuncDef';
for i=1:nargin
    line = varargin{i};
    if(strncmp(line,'-bcType=',7))
        bcType=sscanf(line,'-bcType=%i');
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
    elseif(strncmp(line,'-relaxFactor=',13))
        relaxFactor=sscanf(line,'-relaxFactor=%e');    
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

getOptions=@(solver,tol,bcType,nx,ny,maxIter,resultsName,counter,branch,funcDefFile,ThermalLoading)...
        {'-case=coupledSystem',...
         '-nonlinear',...
         '-saveIC',...
         '-noplot',...
         sprintf('-solver=%s',solver),...
         sprintf('-tol=%e',tol),...
         sprintf('-bcType=%i',bcType),...
         sprintf('-nx=%i',nx),...
         sprintf('-ny=%i',ny),...
         sprintf('-maxIter=%i',maxIter),...
         sprintf('-f=%s%i%c',resultsName,counter,branch),...
         sprintf('-funcDefFile=%s',funcDefFile),...
         sprintf('-xi=%e',ThermalLoading)
         };


% 1st branch (top): from xiMax -> xiMin with step -dxi
counter=0;
ThermalLoading=xiMax-dxi*counter;
exitflag=1;
b=1;
while(ThermalLoading>xiMin && exitflag>0)
    ThermalLoading=xiMax-dxi*counter;
    counter=counter+1;
    options=getOptions(solver1,tol,bcType,nx,ny,maxIter,resultsName,counter,branches(b),funcDefFile,ThermalLoading);
    if(counter>2) % use two previous results as initial guess
        options=[options,...
            {sprintf('-readICResult1=%s',sprintf('%s%i%c',resultsName,counter-1,branches(b))),...
            sprintf('-readICResult2=%s',sprintf('%s%i%c',resultsName,counter-2,branches(b)))}...
            ];
    end
    printCmd(options);
    exitflag = runShell(options{:});
    if(exitflag<=0)
        % this run is not valid. Do not count it and rm the results dir
        resultDir=sprintf('%s%i%c',resultsName,counter,branches(b));
        fprintf('%sThis run is not valid. Do not count it and rm the results dir: %s\n',infoPrefix,resultDir)
        rmdir(resultDir,'s');
        counter=counter-1; 
    end
end
%xiMin=ThermalLoading; % replace xiMin to where the top branch stopped

nTop=counter;
xiTop=ThermalLoading;


% 2nd branch (middle): from xiMin -> xiMax with step dxi
counter=0;
ThermalLoading=xiMin+dxi*counter;
exitflag=1;
b=2;
while(ThermalLoading<xiMax && exitflag>0)
    ThermalLoading=xiMin+dxi*counter;
    counter=counter+1;
    options=getOptions(solver2,tol,bcType,nx,ny,maxIter,resultsName,counter,branches(b),funcDefFile,ThermalLoading);
    if(counter>2) % use two previous results as initial guess
        options=[options,...
            {sprintf('-readICResult1=%s',sprintf('%s%i%c',resultsName,counter-1,branches(b))),...
            sprintf('-readICResult2=%s',sprintf('%s%i%c',resultsName,counter-2,branches(b)))}...
            ];
    end
    printCmd(options);
    exitflag = runShell(options{:});
    if(exitflag<=0)
        % this run is not valid. Do not count it and rm the results dir
        resultDir=sprintf('%s%i%c',resultsName,counter,branches(b));
        fprintf('%sThis run is not valid. Do not count it and rm the results dir: %s\n',infoPrefix,resultDir)
        rmdir(resultDir,'s');
        counter=counter-1; 
    end

end
nMiddle=counter;
xiMiddle=ThermalLoading;


% 3rd branch (lower): from xiMin -> xiMax with step dxi
load(sprintf('%s%i%c/results.mat',resultsName,nTop,branches(1)),'x','xi');
x0=-x; % use the negative results of the top branch as the initial guess for the lower branch
dlmwrite('bifurIC3.dat',x0);

xiMin=xi; % xiMin is where the top branch stopped
counter=0;
ThermalLoading=xiMin+dxi*counter;
exitflag=1;
b=3;

while(ThermalLoading<xiMax && exitflag>0)
    ThermalLoading=xiMin+dxi*counter;
    counter=counter+1;
    options=getOptions(solver3,tol,bcType,nx,ny,maxIter,resultsName,counter,branches(b),funcDefFile,ThermalLoading);
    if(counter>2) % use two previous results as initial guess
        options=[options,...
            {sprintf('-readICResult1=%s',sprintf('%s%i%c',resultsName,counter-1,branches(b))),...
            sprintf('-readICResult2=%s',sprintf('%s%i%c',resultsName,counter-2,branches(b)))}...
            ];
    else
        options=[options,{'-readICFile=bifurIC3'}]; % use ic saved in file
    end
    printCmd(options);
    exitflag = runShell(options{:});
    if(exitflag<=0)
        % this run is not valid. Do not count it and rm the results dir
        resultDir=sprintf('%s%i%c',resultsName,counter,branches(b));
        fprintf('%sThis run is not valid. Do not count it and rm the results dir: %s\n',infoPrefix,resultDir)
        rmdir(resultDir,'s');
        counter=counter-1; 
    end
end
nLower=counter;
xiLower=ThermalLoading;

fprintf('%sTop branch stopped at xi=%e\n',infoPrefix,xiTop);
fprintf('%sNumber results=%i\n',infoPrefix,nTop);
fprintf('%sMiddle branch stopped at xi=%e\n',infoPrefix,xiMiddle);
fprintf('%sNumber results=%i\n',infoPrefix,nMiddle);
fprintf('%sLower branch stopped at xi=%e\n',infoPrefix,xiLower);
fprintf('%sNumber results=%i\n',infoPrefix,nLower);


color='b';
showIC=true;
distBranch=true;
plotBifurcation(resultsName,color,distBranch,showIC);


end



























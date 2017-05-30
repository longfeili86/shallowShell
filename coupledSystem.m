function exitflag=coupledSystem(parameters)
% given parameters, this function solves the coupled system:
%   pde:  \nabla^4 phi = -1/2 L[w,w]*lambda - L[w0,w] -f.phi
%         \nabla^4 w   = L[w,phi]*lambda + L[w0,phi] + f.w
%   (Note: for linear case: lambda=0, for nonlinear case: lambda=1)
%   bcTypes: 0 periodic; 1 simply supported; 2 clamped edge; 3 free edge; 4 CS; 5 CF
%         0: periodic bc  to do FINISH ME ......
%         1: w=0, d^2wdn^2=0, phi=0, d^2phidn^2=0,
%         2: w=0, dwdn=0,  phi=0, dphidn=0    
%         3: d^2wdn^2+nu*d^2ds^2=0, d^3wdn^3+(2-nu)d^3wdnds^2=0, phi=0, dphidn=0 
%         4: mixed CS to do FINISH ME ......
%         5: mixed CF to do FINISH ME ......
%   output:
%           exitflag>0  success
%           exitflag=0  failed interations greater than maxIters
%           exitflag<0  other failure
%
% --Longfei Li

infoPrefix = '--coupledSystem--: '; % all info displayed by this function includes this prefix


% parse parameters
resultsDir=parameters.resultsDir;

xa=parameters.xa;xb=parameters.xb;ya=parameters.ya;yb=parameters.yb;
domain=[xa,xb,ya,yb]; % rectangle [xa,xb,ya,yb]
nu=parameters.nu; % physical parameters (poisson ratio)
bcType=parameters.bcType;% bcType: 2 clamped edge, 3 free edge
nx=parameters.nx; % number of grid points in x direction
ny=parameters.ny; % number of grid points in y direction

isPlot=parameters.isPlot;
savePlot=parameters.savePlot;
saveIC=parameters.saveIC;
useLU=parameters.useLU;

isLinear=parameters.isLinear;
solver=parameters.solver;

relaxFactor=parameters.relaxFactor;
implicitFactor=parameters.implicitFactor;


maxIter=parameters.maxIter;
tol=parameters.tol;

funcDefFile=parameters.funcDefFile;
knownExactSolution=parameters.knownExactSolution;

readICFile=parameters.readICFile; % if non-empty, read IC from this file
saveICFile=parameters.saveICFile; % save the w solution into an ICFile

% if non empty, obtain initial guess from results saved in those dirs
readICResult1=parameters.readICResult1;
readICResult2=parameters.readICResult2;

xi=parameters.xi; %thermal loading

usePAC=parameters.usePAC;

% preprocess

% define grid
myGrid = buildGrid(domain,nx,ny);
Xvec = myGrid.XX(:);%column vector
Yvec = myGrid.YY(:);%column vector
n=length(Xvec); % number of nodes




% define diff matrix
hx = myGrid.hx;
hy = myGrid.hy;
mtx = getDiffMatrix(nx,ny,hx,hy);

% define index
Index=getIndex(nx,ny);

% define a transition function for mixed bc types
% and pass it to parameters. Here we assume the clamped region are the same
% on both top and bottom regions
if(bcType==4 || bcType==5)
    xc=parameters.xc;
    rc=parameters.rc;
    parameters.omega=smoothTrans(Xvec(Index.BoundaryT),xc,rc);
end


% define given functions
isReadIC=(~strcmp(readICResult1,'') || ~strcmp(readICFile,'')); % readIC=true if readICFile is non-empty or if readICResult1 is non-empty
fprintf('%sGetting definitions of all the given functions from file: %s.m\n',infoPrefix,funcDefFile);
run(funcDefFile); % f.w(x,y), f.phi(x,y) and w0(x,y) are defined here
Fphi=f.phi(Xvec,Yvec);
Fw=f.w(Xvec,Yvec);
W0=w0(Xvec,Yvec);
xOld=0.;dx=0.; % this is needed by PAC. Set zero otherwise
if(~usePAC)
    isICFuncDefined=exist('wi','var'); % isICFuncDefined=true if a function_handle wi=@(x,y) is defined in funcDef
    if( isICFuncDefined && ~isReadIC) % get ic from the defined function_handle
        fprintf('%sComputing initial guess from function wi(x,y) defined in %s.m\n',infoPrefix,funcDefFile);   
        Wi = wi(Xvec,Yvec); % evaluate initial guess at nodes. 
    elseif(~isICFuncDefined && isReadIC) % get ic from data file
        if( ~strcmp(readICFile,''))
            ICFileName=sprintf('%s.dat',readICFile);
            fprintf('%sReading in initial guess from data file: %s\n',infoPrefix,ICFileName);      
            x0=dlmread(ICFileName); % the saved ic for all solutions now
        elseif(~strcmp(readICResult1,''))
            x0=readICfromSavedResults(xi,readICResult1,readICResult2);  
        else
            error('something wrong with readIC\n');
        end
        assert((length(x0)==2*n || length(x0)==2*n+3),'Error: size of the read-in initial guess does not match the grid');
    elseif(isICFuncDefined && isReadIC)
        error('Both an IC data file and an IC function are given. I do not know which one to use. Specify only one or none to use W0 by default');
    else
        fprintf('%sUsing w0 as initial guess\n',infoPrefix);           
        Wi=W0;  % if no initial guess is specified, use W0 as initial guess
    end
else
    fprintf('%sUse Pseduo-Arclength Continuation (PAC) Method\n',infoPrefix);
    assert(~strcmp(readICResult1,'') && ~strcmp(readICResult2,'')) % we must provide two previous solutions
    load(sprintf('%s/results.mat',readICResult1),'x','xi');
    xOld=[x;xi];
    load(sprintf('%s/results.mat',readICResult2),'x','xi');
    xOlder=[x;xi];
    dsOld=sqrt(sum((xOld-xOlder).^2));
    fprintf('%scomputed ds=%e\n',infoPrefix,dsOld);
    dx=(xOld-xOlder)/dsOld;
    dx=dx/norm(dx); % normalize it
    x0=xOld+dx*parameters.ds;
end
%RHSs: bc are already implemented inside of getRHS functions
RHSphi=@(w,phi) getRHS_phiEqn(w,phi,Fphi,W0,mtx,parameters,Index);
RHSw=@(w,phi)   getRHS_wEqn(w,phi,Fw,W0,mtx,parameters,Index);

% initial guess for w is known, we obtain initial guess for the whole 
% system x0 by using 1 step of picard:
if(exist('Wi','var')) 
    Aphi = getMTX_phiEqn(Index,mtx,parameters);
    x0=getInitialGuess(n,Wi,Aphi,RHSphi,parameters);
end

% additional RHS for augemented system (free bc)
R =zeros(3,1); 
if(bcType==3)
    if(knownExactSolution)
        % if exact solution is known, then the additional RHS should be given
        % by exact solutions
        addRHS = 0.*Xvec;
        addRHS(Index.UsedPoints)=exact.w(Xvec(Index.UsedPoints),Yvec(Index.UsedPoints));
        Q=getKernalOfSingularMatrix(myGrid,Index);
        R=Q'*addRHS; 
    end
    fprintf('%sFree BC additional rhs: r1=%f;r2=%f;r3=%f\n',infoPrefix,R(1),R(2),R(3));
end




% print some information before solve
sysInfo='nonlinear';
if(isLinear)
    sysInfo='linear';
end
fprintf('%sThe coupled system is %s\n',infoPrefix,sysInfo);
fprintf('%sbcType:  %i\n',infoPrefix,bcType);
fprintf('%sUsing solver: %s\n',infoPrefix,solver);
if(strcmp(solver,'imPicard'))
    fprintf('%simplicitFactor = %f\n',infoPrefix,implicitFactor);  
    assert(implicitFactor>=0 && implicitFactor<=1.,'implicitFactor has to be:  0<=implicitFactor<=1.');
elseif(strcmp(solver,'newton'))
   fprintf('%srelaxFactor = %f\n',infoPrefix,relaxFactor); 
   assert(relaxFactor>0 && relaxFactor<=1.,'relaxFactor has to be: 0<relaxFactor<=1.');
end



% solve the coupled system
if(strcmp(solver,'fsolve'))
    problem=setupFSolveProblem(myGrid,Index,mtx,RHSphi,RHSw,R,W0,x0,n,parameters,xOld,dx);   
    [x,fval,exitflag,output] = fsolve(problem);
elseif(strcmp(solver,'newton'))
    [x,fval,exitflag,output]=newtonSolve(n,x0,W0,Index,mtx,parameters,myGrid,RHSphi,RHSw,R,xOld,dx); % newton solve
elseif(strcmp(solver,'imPicard') || strcmp(solver,'exPicard'))
    [x,fval,exitflag,output]=picardSolve(n,x0,W0,Index,mtx,parameters,myGrid,RHSphi,RHSw,R,xOld,dx); % picard solve
else
    fprintf('%sError unknown solver: %s\n',infoPrefix,solver);
    return
end


% parse solutions
Nphi=1:n; % phi solutions
Nw=n+1:2*n;  % w solutions
Nlambda=2*n+1:2*n+3; % lambda solutions

PHI =x(Nphi); 
W =  x(Nw);
if(bcType==3)
    lambda = x(Nlambda); 
     fprintf('%sFree BC additional variables: lambda1=%e, lambda2=%e, lambda3=%e\n',...
        infoPrefix,lambda(1),lambda(2),lambda(3));
end

if(usePAC)
    xi=x(end);  % the last solution is xi
    x=x(1:end-1);
    fprintf('%sPAC solution xi=%e\n',infoPrefix,xi);
end

% check the solutions at unused points: must be zero
if(false) % everything looks good. No need to check this now
    for i=1:length(Index.UnusedGhostCorners)
        fprintf('%sw Solution at unused point %i: %e\n',infoPrefix,i,W(Index.UnusedGhostCorners(i)));
        fprintf('%sphi Solution at unused point %i: %e\n',infoPrefix,i,PHI(Index.UnusedGhostCorners(i)));
    end
end


% write summary
writeSummary(fval,exitflag,output);

% save solution into an IC file, so that other runs can read as its initial guess
if(~strcmp(saveICFile,''))
    saveICFileName=sprintf('%s.dat',saveICFile);
    fprintf('%sSave solution into an IC file: %s\n',infoPrefix,saveICFileName);
    dlmwrite(saveICFileName,x);
end

% postprocess results
Xplot = reshape(Xvec(Index.interiorBoundary),ny,nx);
Yplot = reshape(Yvec(Index.interiorBoundary),ny,nx);
W0plot = reshape(W0(Index.interiorBoundary),ny,nx);
Wplot = reshape(W(Index.interiorBoundary),ny,nx);
PHIplot = reshape(PHI(Index.interiorBoundary),ny,nx);
Fwplot=reshape(Fw(Index.interiorBoundary),ny,nx);
Fphiplot=reshape(Fphi(Index.interiorBoundary),ny,nx);

% Longfei 20170508: save x into results. So we can load x and use it as
% initial guess for other runs
save(sprintf('%s/results.mat',resultsDir),'Xplot','Yplot','Wplot','W0plot','PHIplot','Fwplot','Fphiplot','x','xi');
if(saveIC)
    PHIi=x0(Nphi);
    Wi=x0(Nw);
    PHIiplot=reshape(PHIi(Index.interiorBoundary),ny,nx); 
    Wiplot=reshape(Wi(Index.interiorBoundary),ny,nx);
    save(sprintf('%s/results.mat',resultsDir),'Wiplot','PHIiplot','-append');
end
if(knownExactSolution)
    WerrPlot=exact.w(Xplot,Yplot)-Wplot;
    PHIerrPlot=exact.phi(Xplot,Yplot)-PHIplot;
    save(sprintf('%s/results.mat',resultsDir),'WerrPlot','PHIerrPlot','-append');
end


if (isPlot)
    figure
    mySurf(Xplot,Yplot,W0plot,'$w_0$');
    if(savePlot)
        printPlot('w0',resultsDir);
    end    
    
    figure
    mySurf(Xplot,Yplot,Wplot,'$w$');
    if(savePlot)
        printPlot('wSolution',resultsDir);
    end
    
    figure
    mySurf(Xplot,Yplot,PHIplot,'$\phi$');
    if(savePlot)
        printPlot('phiSolution',resultsDir);
    end   
    
    figure
    mySurf(Xplot,Yplot,Fwplot,'$f_w$');
    if(savePlot)
        printPlot('wForcing',resultsDir);
    end
    
    figure
    mySurf(Xplot,Yplot,Fphiplot,'$f_{\phi}$');
    if(savePlot)
        printPlot('phiForcing',resultsDir);
    end
    
    % plot error if exact solution is known
    if(knownExactSolution)
        figure 
        mySurf(Xplot,Yplot,WerrPlot,'$E(w)$');
        if(savePlot)
            printPlot('wError',resultsDir);
        end
        figure 
        mySurf(Xplot,Yplot,PHIerrPlot,'$E(\phi)$');
        if(savePlot)
            printPlot('phiError',resultsDir);
        end
    end
end



end
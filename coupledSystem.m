function coupledSystem(parameters)
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
useLU=parameters.useLU;

isLinear=parameters.isLinear;
solver=parameters.solver;

maxIter=parameters.maxIter;
tol=parameters.tol;

funcDefFile=parameters.funcDefFile;
knownExactSolution=parameters.knownExactSolution;




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

% define given functions
fprintf('%sGetting definitions of all the given functions from file: %s.m\n',infoPrefix,funcDefFile);
run(funcDefFile); % f.w(x,y), f.phi(x,y) and w0(x,y) are defined here
Fphi=f.phi(Xvec,Yvec);
Fw=f.w(Xvec,Yvec);
W0=w0(Xvec,Yvec);
PHI0=0.*W0; % store initial guess for phi

% print some information before solve
sysInfo='nonlinear';
if(isLinear)
    sysInfo='linear';
end
fprintf('%sThe coupled system is %s\n',infoPrefix,sysInfo);
fprintf('%sbcType:  %i\n',infoPrefix,bcType);
fprintf('%sUsing solver: %s\n',infoPrefix,solver);




%RHSs: bc are already implemented inside of getRHS functions
RHSphi=@(w,phi) getRHS_phiEqn(w,phi,Fphi,W0,mtx,parameters,Index);
RHSw=@(w,phi)   getRHS_wEqn(w,phi,Fw,W0,mtx,parameters,Index);
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

% solve the coupled system
if(strcmp(solver,'fsolve'))
    problem=setupFSolveProblem(myGrid,Index,mtx,RHSphi,RHSw,R,W0,PHI0,n,parameters);   
    [x,fval,exitflag,output] = fsolve(problem);
elseif(strcmp(solver,'newton'))
    x=newtonSolve(n,PHI0,W0,Index,mtx,parameters,myGrid,RHSphi,RHSw,R); % newton solve
elseif(strcmp(solver,'imPicard') || strcmp(solver,'exPicard'))
    x=picardSolve(n,PHI0,W0,Index,mtx,parameters,myGrid,RHSphi,RHSw,R); % picard solve
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

% check the solutions at unused points: must be zero
if(false) % everything looks good. No need to check this now
    for i=1:length(Index.UnusedGhostCorners)
        fprintf('%s w Solution at unused point %i: %e\n',infoPrefix,i,W(Index.UnusedGhostCorners(i)));
        fprintf('%s phi Solution at unused point %i: %e\n',infoPrefix,i,PHI(Index.UnusedGhostCorners(i)));
    end
end


% postprocess results
Xplot = reshape(Xvec(Index.interiorBoundary),ny,nx);
Yplot = reshape(Yvec(Index.interiorBoundary),ny,nx);
Wplot = reshape(W(Index.interiorBoundary),ny,nx);
PHIplot = reshape(PHI(Index.interiorBoundary),ny,nx);
Fwplot=reshape(Fw(Index.interiorBoundary),ny,nx);
Fphiplot=reshape(Fphi(Index.interiorBoundary),ny,nx);

save(sprintf('%s/results.mat',resultsDir),'Xplot','Yplot','Wplot','PHIplot','Fwplot','Fphiplot');
if(knownExactSolution)
    WerrPlot=exact.w(Xplot,Yplot)-Wplot;
    PHIerrPlot=exact.phi(Xplot,Yplot)-PHIplot;
    save(sprintf('%s/results.mat',resultsDir),'WerrPlot','PHIerrPlot','-append');
end


if (isPlot)
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
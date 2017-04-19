function biharmonic(parameters)
% given parameters, this function solves:
%   pde: \nabla^4 w = rhs
%   bcTypes: 0 periodic; 1 simply supported; 2 clamped edge; 3 free edge; 4 CS; 5 CF
%         0: periodic bc  to do FINISH ME ......
%         1: w=0, d^2wdn^2=0
%         2: w=0, dwdn=0
%         3: d^2wdn^2+nu*d^2ds^2=0, d^3wdn^3+(2-nu)d^3wdnds^2=0
%         4: mixed CS to do FINISH ME ......
%         5: mixed CF to do FINISH ME ......
% --Longfei Li

infoPrefix = '--biharmonic--: '; % all info displayed by this function includes this prefix


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

rhsFile=parameters.rhsFile;
knownExactSolution=parameters.knownExactSolution;

% preprocess

% define grid
myGrid = buildGrid(domain,nx,ny);
Xvec = myGrid.XX(:);%column vector
Yvec = myGrid.YY(:);%column vector

% define diff matrix
hx = myGrid.hx;
hy = myGrid.hy;
mtx = getDiffMatrix(nx,ny,hx,hy);
A = mtx.BiDh;

% define index
Index=getIndex(nx,ny);

% define rhs
fprintf('%sGetting rhs definition from file: %s.m\n',infoPrefix,rhsFile);
run(rhsFile);
RHS = rhs.w(Xvec,Yvec);%vectorized rhs


% solve

% assign bc
A = assignBoundaryConditionsCoefficient(A,Index,mtx,parameters);
RHS = assignBoundaryConditionsRHS(RHS,Index,parameters);
% we need cornor condition for some bc
[A,RHS]=assignCornerConditions(A,RHS,Index,mtx,bcType);



% We might need condition number, eig values ect. of the matrix.
% There are unused ghost points, we need to remove those points from the A.
Aused=A(Index.UsedPoints,Index.UsedPoints);
RHSused=RHS(Index.UsedPoints);

% for Free bc, the system is singular. We use lagrange multiplier method to
% remove sigularity:
if(bcType==3)
    addRHS=0.*RHS;
    if(knownExactSolution)
        addRHS=exact(Xvec(Index.UsedPoints),Yvec(Index.UsedPoints));
    end
    [Aused,RHSused]=removeSingularity(Aused,RHSused,myGrid,Index,addRHS);
end

% we can use LU decomposition of Aused. This is slower for a single
% biharmonic solve, but for iterative methods, we can reuse L, U to 
% to solve the biharmonic system faster than simply backslash Aused every time.
if(useLU)
    fprintf('%suse LU decomposition\n',infoPrefix);
    %tic;
    [L,U]=lu(Aused,0.);
    %toc;
    %tic;
    y = L\RHSused;
    x = U\y;
    %toc;
else
    %tic;
    x=Aused\RHSused;
    %toc;
end
W = 0.*RHS; % zero out stuff to store solution
if(bcType==3)
    W(Index.UsedPoints) = x(1:end-3);
    lambda=x(end-2:end);
    fprintf('%sFree BC additional variables: lambda1=%e, lambda2=%e, lambda3=%e\n',...
        infoPrefix,lambda(1),lambda(2),lambda(3));
else
    W(Index.UsedPoints) = x;
end


% postprocess results
Xplot = reshape(Xvec(Index.interiorBoundary),ny,nx);
Yplot = reshape(Yvec(Index.interiorBoundary),ny,nx);
Wplot = reshape(W(Index.interiorBoundary),ny,nx);
RHSplot=reshape(RHS(Index.interiorBoundary),ny,nx);
save(sprintf('%s/results.mat',resultsDir),'Xplot','Yplot','Wplot','RHSplot');
if(knownExactSolution)
    errPlot=exact(Xplot,Yplot)-Wplot;
    save(sprintf('%s/results.mat',resultsDir),'errPlot','-append');
end


if (isPlot)
    figure
    mySurf(Xplot,Yplot,Wplot,'Solution');
    if(savePlot)
        printPlot('solution',resultsDir);
    end
    
    figure
    mySurf(Xplot,Yplot,RHSplot,'RHS');
    if(savePlot)
        printPlot('RHS',resultsDir);
    end
    % plot error if exact solution is known
    if(knownExactSolution)
        figure 
        mySurf(Xplot,Yplot,errPlot,'Error');
        if(savePlot)
            printPlot('error',resultsDir);
        end
    end
end


end
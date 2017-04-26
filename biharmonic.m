function biharmonic(parameters)
% given parameters, this function solves the harmonic equation:
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

funcDefFile=parameters.funcDefFile;
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

% define index
Index=getIndex(nx,ny);

% define given functions
fprintf('%sGetting definitions of all the given functions from file: %s.m\n',infoPrefix,funcDefFile);
run(funcDefFile);
Fw = f.w(Xvec,Yvec); %vectorized forcing


% setup equations
A = mtx.BiDh;
RHS=Fw;

% assign bc
A = assignBoundaryConditionsCoefficient(A,Index,mtx,parameters);
RHS = assignBoundaryConditionsRHS(RHS,Index,parameters);

% 20170424: removed assignCornerConditions.
% assignBoundaryConditionsCoefficient and assignBoundaryConditionsRHS takes
% care of the corner conditions.
% we need cornor condition for some bc
% [A,RHS]=assignCornerConditions(A,RHS,Index,mtx,bcType); 


% 20170425: removed Aused, and RHSused. We treat the equations at unused
% points as Identity and the RHS as 0. These are done in
% assignBoundaryConditions..
% We might need condition number, eig values ect. of the matrix.
% There are unused ghost points, we need to remove those points from the A.
%Aused=A(Index.UsedPoints,Index.UsedPoints);
%RHSused=RHS(Index.UsedPoints);

% for Free bc, the system is singular. We use lagrange multiplier method to
% remove sigularity:
if(bcType==3)
    % A is the augmented matrix and Q is the kernal
    [A,Q]=removeMatrixSingularity(A,myGrid,Index); 
    R=zeros(3,1); % additional rhs for the augemented system
    if(knownExactSolution)
        addRHS = 0.*RHS;
        addRHS(Index.UsedPoints)=exact.w(Xvec(Index.UsedPoints),Yvec(Index.UsedPoints));
        R=Q'*addRHS;
    end
    fprintf('%sFree BC additional rhs: r1=%f;r2=%f;r3=%f\n',infoPrefix,R(1),R(2),R(3));
    RHS=[RHS;R]; % augmented RHS
end  

% we can use LU decomposition of Aused. This is slower for a single
% biharmonic solve, but for iterative methods, we can reuse L, U to 
% to solve the biharmonic system faster than simply backslash Aused every time.
if(useLU)
    fprintf('%suse LU decomposition\n',infoPrefix);
    %tic;
    [L,U]=lu(A,0.);
    %toc;
    %tic;
    y = L\RHS;
    x = U\y;
    %toc;
else
    %tic;
    x=A\RHS;
    %toc;
end
n=length(Xvec); % number of nodes
W = x(1:n); 
if(bcType==3)
    lambda=x(n+1:n+3);
    fprintf('%sFree BC additional variables: lambda1=%e, lambda2=%e, lambda3=%e\n',...
        infoPrefix,lambda(1),lambda(2),lambda(3));
end

% check the solutions at unused points
for i=1:length(Index.UnusedGhostCorners)
    fprintf('%sSolution at unused point %i: %e\n',infoPrefix,i,W(Index.UnusedGhostCorners(i)));
end



% postprocess results
Xplot = reshape(Xvec(Index.interiorBoundary),ny,nx);
Yplot = reshape(Yvec(Index.interiorBoundary),ny,nx);
Wplot = reshape(W(Index.interiorBoundary),ny,nx);
Fwplot=reshape(Fw(Index.interiorBoundary),ny,nx);
save(sprintf('%s/results.mat',resultsDir),'Xplot','Yplot','Wplot','Fwplot');
if(knownExactSolution)
    WerrPlot=exact.w(Xplot,Yplot)-Wplot;
    save(sprintf('%s/results.mat',resultsDir),'WerrPlot','-append');
end


if (isPlot)
    figure
    mySurf(Xplot,Yplot,Wplot,'$w$');
    if(savePlot)
        printPlot('wSolution',resultsDir);
    end
    
    figure
    mySurf(Xplot,Yplot,Fwplot,'$f_w$');
    if(savePlot)
        printPlot('wForcing',resultsDir);
    end
    % plot error if exact solution is known
    if(knownExactSolution)
        figure 
        mySurf(Xplot,Yplot,WerrPlot,'$E(w)$');
        if(savePlot)
            printPlot('wError',resultsDir);
        end
    end
end


end